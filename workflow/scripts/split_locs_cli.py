import os
import sys
import argparse
from typing import TextIO


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("ont_sunk_pos", help="ONT mapped sunk pos.")
  ap.add_argument("asm_sunk_pos", help="Assembly mapped sunk pos.")
  ap.add_argument("outdir", help="Output dir.")
  ap.add_argument("haplotype", help="Haplotype.")
  args = ap.parse_args()

  outdir = args.outdir
  haplotype = args.haplotype

  os.makedirs(outdir, exist_ok=True)

  # Stream instead of reading into memory ffs
  all_contigs = set()
  ont_fhs: dict[str, TextIO] = {}
  asm_fhs: dict[str, TextIO] = {}

  print(f"Reading ont sunks in {args.ont_sunk_pos}", file=sys.stderr)
  with open(args.ont_sunk_pos, "rt") as ont_fh:
    for line in ont_fh:
      read, read_pos, contig, contig_start, contig_stop = line.strip().split("\t")

      if contig in ont_fhs:
        fh = ont_fhs[contig]
      else:
        fsafe_contig = contig.replace("#", "_").replace("|", "_")
        fname = os.path.join(outdir, f"{fsafe_contig}_{haplotype}.sunkpos")
        fh = open(fname, "wt")
        ont_fhs[contig] = fh

      # Write to file.
      fh.write(line)

      # Store contig
      all_contigs.add(contig)

  for ont_fh in ont_fhs.values():
    ont_fh.close()

  print(f"Finished splitting ont sunks in {args.ont_sunk_pos}", file=sys.stderr)
  
  print(f"Reading asm sunks in {args.asm_sunk_pos}", file=sys.stderr)
  with open(args.asm_sunk_pos, "rt") as asm_fh:
    for line in asm_fh:
      contig, pos_1, sunk, pos_2 = line.strip().split("\t")
      if not contig in all_contigs:
        continue

      if contig in asm_fhs:
        fh = asm_fhs[contig]
      else:
        fsafe_contig = contig.replace("#", "_").replace("|", "_")
        fname = os.path.join(outdir, f"{fsafe_contig}_{haplotype}.loc")
        fh = open(fname, "wt")
        asm_fhs[contig] = fh
      
      # Write to file.
      fh.write(line)

  for asm_fh in asm_fhs.values():
    asm_fh.close()


if __name__ == "__main__":
  raise SystemExit(main())
