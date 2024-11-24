import sys
import argparse
import pandas as pd
from sympy.solvers import solveset
from sympy import Symbol, Reals
import pyranges as pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai")
    ap.add_argument("--bed")
    ap.add_argument("--locs")
    ap.add_argument("--rlen")
    ap.add_argument("--kmer_len")
    ap.add_argument("--outfile")

    args = ap.parse_args()

    # genome_size=3.055e6
    fai = pd.read_csv(
        args.fai,
        sep="\t",
        header=None,
        names=["contig", "length", "cumlen", "3", "4"],
    )
    genome_size = fai["length"].sum() / 1000

    df = pd.read_csv(
        args.bed, header=None, names="Chromosome Start End".split(), sep="\t"
    )
    df["type"] = "gap"
    df3 = df
    df3.reset_index(inplace=True)
    breaks = pr.PyRanges(df3)

    # breaks

    kmermerge = pd.read_csv(
        args.locs, sep="\t", header=None, names=["chrom", "loc", "kmer", "ID"]
    )  # this could be split too
    kmermerge = kmermerge.drop_duplicates(subset="kmer")
    kmermerge["ID2"] = (
        kmermerge["chrom"].astype(str) + ":" + kmermerge["ID"].astype(str)
    )
    kmer_groups = kmermerge.drop_duplicates(subset="ID2")  #
    kmer_groups["dist"] = kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int)
    kmer_groups["dist"] = kmer_groups["dist"].map(lambda x: max(x, 0))

    lendf = pd.read_csv(
        args.rlen,
        sep="\t",
        header=None,
        names=["rname", "len"],
        dtype={"rname": "string", "len": "uint32"},
    )
    lendf.drop_duplicates(inplace=True)
    lendf.set_index("rname", drop=True, inplace=True)
    lendf.sort_values(by="len", ascending=False, inplace=True)
    lendf["bp"] = lendf["len"].cumsum()
    lendf["kbp"] = lendf["len"] / 1000
    lendf["kbp"] = lendf["kbp"].astype(int)

    kbp_cnt = lendf["kbp"].value_counts(sort=False)

    def cov_prob(intsize, readlen, readcnt, basequal=1):  # everything in kbp ?
        if intsize > readlen:
            return 1
        return (1 - (((readlen - intsize) / genome_size) * (basequal**2))) ** (readcnt)

    x = Symbol("x")
    p = 0.94  # read quality
    q = 1 - p
    r = int(args.kmer_len)  # kmer size

    roots2 = solveset(1 - x + q * ((p) ** (r)) * ((x) ** (r + 1)), x, domain=Reals)
    pinv = 1 / p
    roots2 = [x for x in roots2 if x > 0]
    roots2.sort(key=lambda e: abs(pinv - e))
    assert len(roots2) == 2
    roots3 = roots2[1]

    print("Should be unequal: ", roots3, pinv, file=sys.stderr)

    n = 30  # group size
    qn = ((1 - p * roots3) / (q * (r + 1 - r * roots3))) * (1 / (roots3 ** (n + 1)))
    pn = float(1 - qn)
    # print(qn, pn) # odds perfect 20mer not seen in kmergroup length 22

    covprobs = []
    for intsize in range(1, 3500):
        cum_prob = 1
        for i, x in kbp_cnt.items():
            prob = cov_prob(intsize, i, x, pn)
            cum_prob = cum_prob * prob
        covprobs.append((intsize, 1 - cum_prob))

    covprobsdict = {k: v for k, v in covprobs}
    covprobsdict[0] = 1.0

    for cov, prob in covprobs:
        if prob <= 0.99:
            print(cov, prob)
            break

    max_gaps = []
    max_gaps2 = []
    for x in breaks.df.itertuples():
        ix, ix2, contig, regstart, regend, types = x  # noqa: F841
        xmin = regstart  # noqa: F841
        xmax = regend  # noqa: F841
        kmers = kmer_groups.query(
            "(chrom == @contig) & (ID > (@xmin-2)) & (ID < (@xmax+2))"
        )
        max_gaps.append(kmers["dist"].max())
        max_gaps2.append(
            (
                contig,
                regstart,
                regend,
                kmers["dist"].max(),
                kmers.iloc[kmers["dist"].argmax()].kmer,
            )
        )

    # breaks
    df1 = breaks.df
    df1["max_gap"] = max_gaps
    df1["covprob"] = df1["max_gap"].apply(lambda x: covprobsdict[int(x / 1000)])
    df1.to_csv(args.outfile, index=False, sep="\t")


if __name__ == "__main__":
    main()
