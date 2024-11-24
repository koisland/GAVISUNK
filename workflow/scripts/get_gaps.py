import sys
import pandas as pd
import numpy as np
import argparse
import pyranges as pr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fai")
    parser.add_argument("indir")
    parser.add_argument("outfile_no_data")
    parser.add_argument("outfile_gaps")

    args = parser.parse_args()

    indir = args.indir
    outfile_no_data = args.outfile_no_data
    outfile_gaps = args.outfile_gaps

    fai = pd.read_csv(
        args.fai, sep="\t", header=None, names=["contig", "length", "cumlen", "3", "4"]
    )
    contigs = fai.contig.tolist()
    allbreaks = {}
    outbreaks = []
    bp_breaks = 0
    bp_no_data = 0
    nodata = []
    for x in contigs:
        cont2 = x.replace("#", "_")
        contiglen = fai.query("contig == @x").length.values[0]
        try:
            bedin = pd.read_csv(
                indir + cont2 + ".bed",
                sep="\t",
                header=None,
                names=["contig", "start", "end"],
            )

        except Exception:
            print("No bed for", x, file=sys.stderr)
            nodata.append((x, 0, contiglen))
            bp_no_data += contiglen
            continue
        if len(bedin) == 0:
            nodata.append((x, 0, contiglen))
            bp_no_data += contiglen
            continue

        bed = pr.PyRanges(
            bedin.rename(
                columns={"contig": "Chromosome", "start": "Start", "end": "End"}
            )
        )
        bedin2 = bed.merge().df
        bedin2.columns = ["contig", "start", "end"]
        bedin2["end_prev"] = bedin2["end"].shift(1, fill_value=np.nan)

        breaks = []
        first = True
        for row in bedin2.itertuples():
            if row.end < row.end_prev:
                print("weird", x, file=sys.stderr)
                break
            if first:
                first = False
                # print("confirmed region begins: ",row.start)
            else:
                # print(int(row.end_prev), row.start,)
                breaks.append(
                    (
                        int(row.end_prev),
                        row.start,
                    )
                )
                outbreaks.append(
                    (
                        x,
                        int(row.end_prev),
                        row.start - 1,
                    )
                )
                bp_breaks += row.start - int(row.end_prev)

        #     print(row.end)
        #     print("unconfirmed bp at contig end: ", contiglen - row.end)
        # #     print(breaks)
        allbreaks[x] = breaks

    print("breaks: ", len(outbreaks), bp_breaks, file=sys.stderr)
    print("no data: ", len(nodata), bp_no_data, file=sys.stderr)
    pd.DataFrame(nodata).to_csv(outfile_no_data, header=False, sep="\t", index=False)
    pd.DataFrame(outbreaks).to_csv(outfile_gaps, header=False, sep="\t", index=False)


if __name__ == "__main__":
    main()
