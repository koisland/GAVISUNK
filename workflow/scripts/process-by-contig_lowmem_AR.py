#! /usr/bin/python
import graph_tool as gt
import graph_tool.topology

from scipy.spatial.distance import pdist
import itertools as it

import pandas as pd
import numpy as np

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("SUNKs", help=".loc file of SUNKs aligned to assembly")
    parser.add_argument("sunkpos", help=".sunkpos file of SUNK locations on ONT reads")
    parser.add_argument("rlen", help=".rlen file of ONT read lengths")
    parser.add_argument("badsunks", help="list of SUNKs to exclude")
    parser.add_argument("outputfile", help="output file for intermediate")
    parser.add_argument("outputbed", help="output file for bed")

    parser.add_argument(
        "--minlen",
        help="minimum length filter for reads to visualize",
        default=10000,
        type=int,
    )
    parser.add_argument("--opt_filt", help="apply optional filter", action="store_true")

    args = parser.parse_args()
    # contig = args.contig
    ofile = args.outputfile

    lendf = pd.read_csv(
        args.rlen,
        sep="\t",
        header=None,
        names=["rname", "len"],
        dtype={"rname": "string", "len": "uint32"},
    )
    lendf.set_index("rname", drop=True, inplace=True)
    minlen = args.minlen

    kmermergefile = args.SUNKs
    kmermerge = pd.read_csv(
        kmermergefile, sep="\t", header=None, names=["chrom", "loc", "kmer", "ID"]
    )  # this could be split too
    kmermerge = kmermerge.drop_duplicates(subset="kmer")
    kmermerge["ID2"] = (
        kmermerge["chrom"].astype(str) + ":" + kmermerge["ID"].astype(str)
    )
    sunkposfile = args.sunkpos
    sunkposcat = pd.read_csv(
        sunkposfile,
        sep="\t",
        header=None,
        names=["rname", "pos", "chrom", "start", "ID"],
        dtype={
            "rname": "string",
            "pos": "uint32",
            "chrom": "category",
            "start": "uint32",
            "ID": "uint32",
        },
    )  #

    contig = sunkposcat["chrom"].unique().tolist()[0]

    ## get list of bad sunks
    with open(args.badsunks, "r") as f:
        badsunkin = set(f.read().splitlines())  # noqa: F841
    sunkposcat["ID2"] = (
        sunkposcat["chrom"].astype(str) + ":" + sunkposcat["ID"].astype(str)
    )

    sunkposcat = sunkposcat.query("ID2 not in @badsunkin")

    kmer_groups = kmermerge.drop_duplicates(subset="ID")  #
    pd.options.mode.chained_assignment = None  # default='warn
    kmer_groups["dist"] = (
        kmer_groups.ID.diff().fillna(kmer_groups.ID).astype(int)
    )  # rows)

    # Filter 1: only keep reads with at least two SUNKs within target region
    multisunk = (
        sunkposcat.groupby("rname", as_index=False)
        .agg({"ID": "nunique"})
        .sort_values(by="ID")
        .query("ID > 1")
    )
    if len(multisunk) == 0:
        pd.DataFrame(list([contig])).to_csv(ofile, header=False, sep="\t", index=False)
        exit()

    ### RESTRICT TO REGION
    multisunkset = set(multisunk.rname.tolist())  # noqa: F841
    sunkpos3 = sunkposcat.query("rname in @multisunkset ").sort_values(
        by=["rname", "start"]
    )  # reads with multiple SUNK group matches # & start >= @xmin & start <= @xmax

    sunkpos3b = sunkpos3.drop_duplicates()

    minlen = 10000  # noqa: F841
    minlenset = set(lendf.query("len >= @minlen").index.tolist())  # noqa: F841
    sub = sunkpos3b.query("rname in @minlenset")  # sunkpos3

    # part2 bookmark

    ### keeep pos

    sub["start"] = sub["start"].astype("int64")
    sub["pos"] = sub["pos"].astype("int64")

    outputs = []
    counter = 0
    sub["start"] = sub["start"].astype("int64")
    sub["pos"] = sub["pos"].astype("int64")
    grouped = sub.groupby("rname")
    for rname, g in grouped:
        counter += 1
        startarray = pdist(g[["start"]])
        posarray = pdist(g[["pos"]])
        signarray = pdist(g[["pos"]], lambda u, v: u[0] > v[0])
        idarray = list(it.combinations(g["ID"], 2))
        locarray = list(
            it.combinations(g["pos"], 2)
        )  # keep track of ID-Pos corresponence
        with np.errstate(
            divide="ignore"
        ):  # kmers occasionally seen in multiple locations on same read
            diffarray = posarray / startarray
        mask = np.logical_not((diffarray >= 1.1) + (diffarray <= 0.9))
        if sum(mask) < 1:
            continue

        # only calculate on masked version (distdf2 equiv)
        signarray = signarray.astype(int)
        vals, counts = np.unique(signarray[mask], return_counts=True)
        trueOrient = vals[np.argmax(counts)]
        mask2 = np.logical_and(mask, (signarray == trueOrient))

        stack = np.column_stack([idarray, locarray])[mask2, :]
        sub2 = pd.DataFrame(np.column_stack([idarray, locarray])[mask2, :])
        sub2.columns = ["ID", "ID2", "pos1", "pos2"]
        mergedf = pd.DataFrame(np.row_stack([stack[:, [0, 2]], stack[:, [1, 3]]]))
        mergedf.columns = ["ID", "pos"]

        mergedf = pd.concat(
            [
                sub2[["ID", "pos1"]].rename(columns={"pos1": "pos"}),
                sub2[["ID2", "pos2"]].rename(columns={"ID2": "ID", "pos2": "pos"}),
            ]
        )

        multipos = (
            mergedf.groupby("ID")["pos"].agg(nunique="nunique").query("nunique>1")
        )

        if len(multipos) >= 1:
            badrowlist = []
            for i in multipos.index:
                goodpos = (  # noqa: F841
                    mergedf.query("ID == @i")["pos"].value_counts().index[0]
                )  # any checks that this is good?
                badrowlist += list(
                    sub2.query(
                        "(ID == @i) & not ((pos1 == @goodpos)|(pos2 == @goodpos))"
                    ).index
                )

            sub2b = sub2.drop(badrowlist)
        else:
            sub2b = sub2

        g = gt.Graph(directed=False)
        name = g.add_edge_list(
            sub2b[["ID", "ID2"]].values, hashed=True, hash_type="int64_t"
        )
        g.vp.name = name
        glarge = graph_tool.topology.extract_largest_component(g)

        glarge.vp.name.get_array()
        df = pd.DataFrame(g.vp.name.get_array()[glarge.get_vertices()], columns=["ID"])
        df["rname"] = rname
        outputs.append(df)

    if len(outputs) == 0:
        pd.DataFrame(list([contig])).to_csv(ofile, header=False, sep="\t", index=False)
        exit()
    else:
        outputsdf = pd.concat(outputs)
        outputsdf.to_csv(ofile, header=False, sep="\t", index=False)

    # made outputdf

    ### make contig-wide graph
    idarray = pd.DataFrame(
        outputsdf.groupby("rname").apply(lambda g: list(it.combinations(g["ID"], 2)))
    ).explode([0])
    idarray.rename(
        columns={
            0: "IDs",
        },
        inplace=True,
    )
    idarray[["ID", "ID2"]] = pd.DataFrame(idarray["IDs"].tolist(), index=idarray.index)

    g = gt.Graph(directed=False)
    name = g.add_edge_list(
        idarray[["ID", "ID2"]].values, hashed=True, hash_type="int64_t"
    )
    g.vp.name = name

    c = graph_tool.topology.label_components(g)[0]

    regions = []
    for s in list(pd.unique(c.a)):
        sunks = g.vp.name.get_array()[gt.GraphView(g, vfilt=c.a == s).get_vertices()]
        if len(sunks) <= 2:
            continue
        start = int(sunks.min())
        end = int(sunks.max())
        span = end - start
        regions.append((start, end, span, len(sunks)))

    regdf = pd.DataFrame(
        regions, columns=["start", "end", "span", "sunks"]
    ).sort_values(by="start")

    regdf["contig"] = contig
    regdf[["contig", "start", "end"]].to_csv(
        args.outputbed, header=False, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
