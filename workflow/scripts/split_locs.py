import pandas as pd
import os
import argparse
from pandas.core.groupby.generic import DataFrameGroupBy


def splitlocs(group_df: DataFrameGroupBy[tuple], dirname: str, exten: str):
    for name, contig_group in group_df:
        name = name[0]
        contig_name = name.replace("#", "_")
        output_path = f"{dirname}/{contig_name}.{exten}"
        contig_group.to_csv(output_path, header=None, index=None, sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ont_pos", required=True)
    ap.add_argument("--kmer_loc", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    ONT_pos_df = pd.read_csv(
        args.ont_pos,
        header=None,
        sep="\t",
        names=["read", "read_pos", "contig", "contig_start", "contig_stop"],
    )
    kmer_loc = pd.read_csv(
        args.kmer_loc,
        header=None,
        sep="\t",
        names=["contig", "pos1", "SUNK", "pos2"],
    )
    contig_list = ONT_pos_df["contig"].unique().tolist()
    kmer_loc_subset = kmer_loc[kmer_loc["contig"].isin(contig_list)]
    group_ONT = ONT_pos_df.groupby(["contig"])
    group_loc = kmer_loc_subset.groupby(["contig"])
    splitlocs(group_ONT, outdir, "sunkpos")
    splitlocs(group_loc, outdir, "loc")


if __name__ == "__main__":
    main()
