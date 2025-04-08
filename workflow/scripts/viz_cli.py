import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Rectangle
import matplotlib.ticker as mtick
import argparse
import os
from ncls import NCLS
from concurrent.futures import ProcessPoolExecutor


def plot_contig(contig, args):
    cont2 = contig.replace("#","_")
    print(contig)
    bedreg= bedreg1.query('chrom == @contig')

    try: outputsdf= pd.read_csv(args.rgraph_dir + "/" + cont2 + "_" + args.haplotype + ".tsv" ,sep="\t",names=['ID','rname'],dtype = {'ID':'int','rname':'string'})
    except:
        print("outputsdf not found for", contig)

    outputsdf.set_index("rname",inplace=True,drop=True)

    if len(outputsdf) <2:
        print("NO READS FOR", contig)


    kmermergefile = args.rsplitsunk_dir +"/"+ cont2 + "_" + args.haplotype + ".loc"
    kmermerge = pd.read_csv(kmermergefile,sep="\t",header=None,names=['chrom','loc','kmer','ID'])# this could be split too
    kmermerge = kmermerge.drop_duplicates(subset='kmer')

    sunkposfile = args.rsplitsunk_dir+"/"+ cont2 + "_" + args.haplotype + ".sunkpos"
    sunkposcat = pd.read_csv(sunkposfile,sep="\t",header=None,names=['rname','pos','chrom','start','ID'],dtype={'rname':'string','pos':'int64','chrom':'category','start':'int64','ID':'int64'})

    kmermerge['ID2'] = kmermerge['chrom'].astype(str) +":"+ kmermerge['ID'].astype(str)
    sunkposcat['ID2'] = sunkposcat['chrom'].astype(str) +":"+ sunkposcat['ID'].astype(str)
    print("sunkposcat len: ", len(sunkposcat))
  
    kmer_groups = kmermerge.drop_duplicates(subset='ID') #  
    kmer_groups['dist'] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int) #rows)
    
    multisunk = sunkposcat.groupby('rname',as_index=False).agg({'ID':'nunique'}).sort_values(by='ID').query('ID > 1')
    print("Filter 1 multisunk len: ", len(multisunk))
    
    if len(multisunk) == 0:
      print("No viable reads for this region: "+str(contig)+"_"+str(regstart)+"_"+str(regend))

    
    multisunkset = set(multisunk.rname.tolist())
    
    ### RESTRICT TO REGION
    sunkpos3 = sunkposcat.query('rname in @multisunkset ').sort_values(by=['rname','start']) # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax

    print(len(sunkpos3))
    print(len(pd.unique(sunkpos3.rname)))
    print(len(pd.unique(sunkposcat.rname)))

    sunkpos3.set_index("rname",drop=True,inplace=True)

    for row in bedreg.itertuples():
        contig = row.chrom
        regstart = row.start
        regend = row.end
        xmin = regstart
        xmax = regend
        breakloc = int((xmin + xmax)/2)
        rightset = set(outputsdf.query("ID > @breakloc").index)
        df = outputsdf.query("ID > @xmin & ID < @xmax")

        if len(df) == 0:
            print("NO HITS FOR", row)
            continue

        kmerlist = list(pd.unique(kmer_groups.query("ID > @xmin & ID < @xmax").ID))
        kmerset = set(kmerlist)
        never_seen2 = kmerset - set(df.ID)
        print(len(never_seen2))
        constart = df.ID.min()
        conend = df.ID.max()
        maxes={}
        first=True
        fig,ax = plt.subplots(figsize=(20,9))
        fig.subplots_adjust(right=0.75)
        ax.set_xlim(regstart,regend)
        intervals = []
        poss = []
        print("NUMBER OF RNAMES: ",len(list(pd.unique(df.index))))
        for rname in list(pd.unique(df.index)):
            try: rlen = lendf.loc[rname].len
            except:
                print("no rlen for: ", rname)
                if runmode == 'user_bed':
                    continue
            sub = df.loc[rname]
            
            try: idlist = set(sub.ID)
            except: idlist = set([sub.ID])

            possub = sunkpos3.loc[rname].query("ID in @idlist")
            posfirst = possub.iloc[0]
            poslast = possub.iloc[-1]
            forward = True
            badlocs_trans = []

            if posfirst.pos > poslast.pos: forward=False
            
            if forward:
                posstart = posfirst.start - posfirst.pos
                posend =  poslast.start + (rlen - poslast.pos)  
            else:
                posend = poslast.pos + poslast.start
                posstart =  posfirst.start - (rlen - posfirst.pos)

            start = posstart
            end = posend

            
            if runmode == 'gaps':
                plotgroup=1
                columns = ['Start','End','Plotgroup']
                if rname in rightset: plotgroup=0
                intervals.append((start,end,plotgroup))
            elif runmode == 'user_bed':
                intervals.append((start,end))
                columns = ['Start','End']

            poss.append(list(idlist))

        interval_df = pd.DataFrame(intervals,columns=columns)
        ncls = NCLS(interval_df.Start,interval_df.End,interval_df.index)
        row_asn = dict()
        interval_df.sort_values(by='Start',inplace=True)
        print("INTERVAL_DF len: ", len(interval_df))
        rows_used = set([0])
        if runmode == 'user_bed':
            for i,s,e in interval_df.itertuples():
                overlap = ncls.find_overlap(s,e)
                rows_overlap = set()
                for i2,s2,e2 in overlap:
                    rows_overlap.add( row_asn.get(e2,0))
                rows_avail = rows_used - rows_overlap
                if len(rows_avail) == 0:
                    row_pick = sorted(rows_used)[-1]+1
                    rows_used.add(row_pick)
                    row_asn[i] = row_pick
                else:
                    row_asn[i] = sorted(rows_avail)[0]
            interval_df['row_asn'] = interval_df.index.map(lambda x: row_asn[x])
            patches = []
            ticks = []
            for i,s,e,r in interval_df.itertuples():
                patches.append(Rectangle((s, r), e-s, 0.7))
                tick = [Rectangle((x, r), 1, 0.7)  for x in poss[i]]
                ticks = ticks + tick
        elif runmode == 'gaps':
            for pgi in [0,1]:
                intv_sub = interval_df.query("Plotgroup == @pgi")
                intv_sub['End'] = intv_sub['End'] + 500 # Minimum separation between reads
                if pgi==1:
                    intv_sub = intv_sub.sort_values(by=['End'],ascending=False)
                for i,s,e,pg in intv_sub.itertuples():

                    overlap = ncls.find_overlap(s,e)
                    rows_overlap = set()
                    for i2,s2,e2 in overlap:
                        rows_overlap.add( row_asn.get(e2,0))

                    rows_avail = rows_used - rows_overlap
                    if len(rows_avail) == 0:
                        row_pick = sorted(rows_used)[-1]+1
                        rows_used.add(row_pick)
                        row_asn[i] = row_pick
                    else:
                        row_asn[i] = sorted(rows_avail)[0]
                    maxrow_prev = max(row_asn.values())

            interval_df['row_asn'] = interval_df.index.map(lambda x: row_asn[x])
            patches = []
            ticks = []
            for i,s,e,pg,r in interval_df.itertuples():
                patches.append(Rectangle((s, r), e-s, 0.7))
                tick = [Rectangle((x, r), 1, 0.7)  for x in poss[i]]
                ticks = ticks + tick

        never_seen3 = [x for x in never_seen2 if not ((x<constart) | (x>conend))]

        if args.colorbed:
            patches2=[]
            for r in dm[['chrStart','chrEnd','color','chr']].query('chr==@contig & chrEnd > @xmin & chrStart < @xmax ').itertuples():
                polygon = Polygon([[r[1],-4],[r[1],-5],[r[2],-5],[r[2],-4]], closed=True, color=r[3],alpha=1)
                patches2.append(polygon)
            p = PatchCollection(patches2, alpha=1, match_original=True)
            ax.add_collection(p)
        ax.add_collection(PatchCollection(patches,color='lightgray',zorder=1,aa=True,edgecolor='w',linewidth=0.01))
        ax.add_collection(PatchCollection(ticks,color='k',zorder=7,linewidth=0.5,edgecolor='k',aa=True))
        ax.scatter(x=never_seen3,y=[-2.5]*len(never_seen3),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
        ax.scatter(x=kmerlist,y=[-1.5]*len(kmerlist),color='k',s=55,marker="|",facecolors=None,linewidths=0.4,zorder=2)
        fmt = '{x:,.0f}'
        tick = mtick.StrMethodFormatter(fmt)
        ax.xaxis.set_major_formatter(tick) 
        ax.set_ylabel("ONT Read Depth")
        ax.set_xlabel("Contig coordinate")
        plt.rcParams['svg.fonttype'] = 'none'
        print(outdir + "/" + args.sample +"_" + args.haplotype + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".svg")
        plt.savefig(outdir + "/" + args.sample +"_" + args.haplotype + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".svg",format="svg",pad_inches=0,bbox_inches='tight')
        plt.savefig(outdir + "/" + args.sample +"_" + args.haplotype + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".png",pad_inches=0,bbox_inches='tight',dpi=300)
        plt.savefig(outdir + "/" + args.sample +"_" + args.haplotype + "_" +cont2+"_"+str(regstart)+"_"+str(regend)+".pdf",pad_inches=0,bbox_inches='tight')
        plt.close()
        print("Plotting complete for " +contig+"_"+str(regstart)+"_"+str(regend))

        output_sunks_fh = open(f"{outdir}/{contig}_sunks.txt", "wt")
        for kmer in kmerlist:
            print(kmer, file=output_sunks_fh)
        output_sunks_fh.close()

ap = argparse.ArgumentParser()
ap.add_argument("--regbedfile",required=True, type=str, help="Region bedfile to evaluate.")
ap.add_argument("--runmode", required=True, type=str, choices=["gap", "user_bed"])
ap.add_argument("--rsplitsunk_dir", required=True, type=str, help="Read sunks dir split by contig. Expects {read}_{hap}.loc and {read}_{hap}.sunkpos")
ap.add_argument("--rgraph_dir", required=True, type=str, help="Read graph dir. Expects {read}_{hap}.tsv")
ap.add_argument("--rlen", required=True, type=str, help="Read lengths by haplotype.")
ap.add_argument("--colorbed", default=None, help="Optional color bed file.")
ap.add_argument("--sample", required=True)
ap.add_argument("--haplotype", required=True)
ap.add_argument("--outdir", required=True, help="Output directory.")
ap.add_argument("--processes", type=int, default=4, help="Number of proceses to spawn.")
args = ap.parse_args()

runmode=args.runmode

regbedfile = args.regbedfile
if os.stat(regbedfile).st_size == 0:
    print("empty gap file")
    exit()

bedreg1 = pd.read_csv(regbedfile, delimiter="\t",encoding='utf-8',header=None)
header = ['chrom','start','end']
bedreg1.columns = header + [''] * (len(bedreg1.columns) - len(header))  

lendf = pd.read_csv(args.rlen,sep="\t",header=None,names=['rname','len'],dtype={'rname':'string','len':'uint32'})
lendf.drop_duplicates(inplace=True)
lendf.set_index("rname",drop=True,inplace=True)


bedreg1['chrom'].value_counts()
contigs = list(bedreg1['chrom'].value_counts().index)
print(len(contigs))

if args.colorbed:
    dm = pd.read_csv(args.colorbed, delim_whitespace=True, names=['chr','chrStart','chrEnd','color'], header=None,dtype = {'chr':'string','chrStart':'int','chrEnd':'int','color':'string'})
    dm = dm[dm.chrStart >= 0]
    dm["y"] = 1
    dm["func"] = ""

outdir = args.outdir
os.makedirs(outdir, exist_ok=True)


with ProcessPoolExecutor(max_workers=args.processes) as pool:
    pool.map(
        plot_contig,
        *zip(
            *[
                (
                    contig,
                    args
                )
                for contig in contigs
            ]
        ),    
    )

