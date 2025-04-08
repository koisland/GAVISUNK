use std::{
    fs::File, hash::Hash, io::{BufRead, BufReader, BufWriter, Write}, path::PathBuf
};

use clap::Parser;
use itertools::Itertools;
use indexmap::{IndexMap, IndexSet};
use rustc_hash::FxBuildHasher;


type FxIndexMap<K, V> = IndexMap<K, V, FxBuildHasher>;
type FxIndexSet<T> = IndexSet<T, FxBuildHasher>;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(short = 'i', long, required = true)]
    infile: PathBuf,
    #[arg(short = 's', long, required = true)]
    sunk_db: PathBuf,
    #[arg(short = 'o', long, required = true)]
    outdir: PathBuf,
    #[arg(short = 'k', long, default_value_t = 20)]
    kmer_len: usize,
    #[arg(short = 'd', long, default_value_t = 0.02)]
    sunk_dst: f32,
}

#[derive(Eq, PartialEq)]
struct ID {
    read: Box<str>,
    id: usize,
}
impl Hash for ID {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.read.hash(state);
        self.id.hash(state);
    }
}
fn main() -> eyre::Result<()> {
    let cli = Cli::parse();

    let prefix = cli.infile.file_stem().unwrap();

    let _ = std::fs::create_dir(&cli.outdir);

    let pair_file = cli.outdir.join(prefix).with_extension("pair");
    // let dst_file = cli.outdir.join(prefix).with_extension("dst");
    let interm_file = cli.outdir.join(prefix).with_extension("sunkpos");
    let interm_cnt_file = cli.outdir.join(prefix).with_extension("sunknum");

    let mut pair_writer = BufWriter::new(File::create(pair_file)?);
    // let mut dst_writer = BufWriter::new(File::create(&dst_file)?);
    let mut interm_writer = BufWriter::new(File::create(&interm_file)?);
    let mut interm_cnt_writer = BufWriter::new(File::create(interm_cnt_file)?);

    let sunk_db_reader = BufReader::new(File::open(&cli.sunk_db)?);
    
    eprintln!("Building sunk pool.");
    let mut sunk_info: FxIndexMap<Box<str>, (u16, usize, usize)> = FxIndexMap::default();
    let mut id_ctg_key: FxIndexMap<u16, String> = FxIndexMap::default();
    let mut ctg_id_key: FxIndexMap<String, u16> = FxIndexMap::default();
    let mut ctg_id = 0;
    for line in sunk_db_reader.lines().map_while(Result::ok) {
        let Some((ctg, pos, sunk, id)) = line.trim().split('\t').collect_tuple() else {
            eprintln!("Invalid line. {line}");
            continue;
        };

        let id = id.parse::<usize>().unwrap();
        let pos = pos.parse::<usize>().unwrap();
        let rc_sunk = sunk
            .chars()
            .rev()
            .map(|nt| match nt {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                'N' => 'N',
                _ => panic!("Invalid nt {nt}"),
            })
            .collect::<String>();
        let ctg_id = if let Some(id) = ctg_id_key.get(ctg) {
            *id
        } else {
            ctg_id += 1;
            ctg_id_key.entry(ctg.to_owned()).or_insert(ctg_id);
            id_ctg_key.entry(ctg_id).or_insert(ctg.to_owned());
            ctg_id
        };
        sunk_info.insert(sunk.into(), (ctg_id, pos, id));
        sunk_info.insert(rc_sunk.into(), (ctg_id, pos, id));
    }

    eprintln!("Annotating ONT reads");

    let reads_fa_reader = BufReader::new(File::open(&cli.infile)?);
    let mut fa = noodles::fasta::Reader::new(reads_fa_reader);

    for seq in fa.records().flatten() {
        let read = std::str::from_utf8(seq.name()).unwrap();
        let sequence = std::str::from_utf8(seq.sequence().as_ref())?;
        for (rpos, kmer) in (0..sequence.len() - cli.kmer_len)
            .filter_map(|i| sequence.get(i..i + cli.kmer_len).map(|kmer| (i, kmer)))
        {
            if let Some((ctg, cpos, cid)) = sunk_info
                .get(kmer)
                .map(|(ctg_id, pos, id)| (&id_ctg_key[ctg_id], pos, id))
            {
                writeln!(
                    interm_writer,
                    "{}\t{}\t{}\t{}\t{}",
                    read.to_owned(),
                    rpos,
                    ctg.to_owned(),
                    *cpos,
                    *cid,
                )
                .unwrap();
            }
        }
    }

    let mut seen_ids: FxIndexSet<ID> = FxIndexSet::default();
    let mut read_sunks: FxIndexMap<Box<str>, FxIndexMap<Box<str>, Vec<(usize, usize)>>> =
        FxIndexMap::default();
    let interm_reader = BufReader::new(File::open(&interm_file)?);

    for line in interm_reader.lines().map_while(Result::ok) {
        let Some((read, rpos, ctg, _cst, cpos)) = line.trim().split('\t').collect_tuple() else {
            continue;
        };
        let rpos = rpos.parse::<usize>()?;
        let cpos = cpos.parse::<usize>()?;
        let id = ID {
            read: read.into(),
            id: cpos,
        };
        if seen_ids.contains(&id) {
            continue;
        } else {
            seen_ids.insert(id);
        }

        read_sunks
            .entry(read.into())
            .and_modify(|a| {
                a.entry(ctg.into())
                    .and_modify(|a| a.push((rpos, cpos)))
                    .or_insert_with(|| vec![(rpos, cpos)]);
            })
            .or_insert_with(|| FxIndexMap::from_iter([(ctg.into(), vec![(rpos, cpos)])]));
    }

    eprintln!("Building sunk distance.");
    let mut read_sunk_max: FxIndexMap<&str, FxIndexMap<&str, usize>> = FxIndexMap::default();

    for (read, ctg_pos) in read_sunks.iter() {
        for (ctg, records) in ctg_pos.iter().filter(|(_, rec)| rec.len() > 1) {
            //   2 1 0
            // 0 x x x
            // 1 x x
            // 2 x
            let mut kmer_set = FxIndexSet::default();
            let mut dsts = vec![];
            for (x, (pos_reads, pos_ctg)) in records.iter().enumerate() {
                for (y, (pos_reads_2, pos_ctg_2)) in
                    (0..(records.len() - x - 1)).map(|y| (x + y + 1, records[x + y + 1]))
                {
                    let dst_read = pos_reads_2 as isize - *pos_reads as isize;
                    let dst_ctg = pos_ctg_2 as isize - *pos_ctg as isize;
                    dsts.push((x, y, dst_read, dst_ctg));
                }
            }
            for (x, y, dst_read, dst_ctg) in dsts {
                let dst = dst_ctg.abs_diff(dst_read) as f32 / dst_read as f32;
                // writeln!(
                //     &mut dst_writer,
                //     "{read}\t{ctg}\t{x}\t{y}\t{dst_read}\t{dst_ctg}\t{dst}"
                // )?;
                if dst < cli.sunk_dst {
                    kmer_set.insert(x);
                    kmer_set.insert(y);
                }
            }

            if kmer_set.is_empty() {
                continue;
            }
            read_sunk_max
                .entry(read)
                .and_modify(|n| {
                    n.insert(ctg.as_ref(), kmer_set.len());
                })
                .or_insert_with(|| FxIndexMap::from_iter([(ctg.as_ref(), kmer_set.len())]));
        }
    }
    for (read, sunks) in read_sunk_max.into_iter() {
        // Replace the max sunk.
        let (ctg, cnt) = sunks.iter().max_by(|a, b| a.1.cmp(b.1)).unwrap();
        writeln!(&mut pair_writer, "{read}\t{ctg}\t{cnt}")?;
        for (ctg, cnt) in sunks {
            writeln!(&mut interm_cnt_writer, "{read}\t{ctg}\t{cnt}")?;
        }
    }

    Ok(())
}
