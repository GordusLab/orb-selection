"""Select contamination-flag candidates from the DIAMOND scans of the suspect samples.

Reads one or more DIAMOND tables (outfmt 6: qseqid sseqid pident length qcovhsp scovhsp evalue bitscore;
IDs are '<species>__<gene>'), keeps matches to a species in a different family with protein identity >= 98%
and an alignment >= 100 aa, counts k = the species sharing each suspect gene (itself + distinct partner species),
and keeps the best partner per suspect gene (highest identity, then longest alignment).

Writes in --out-dir:
  candidates.tsv       one row per suspect gene (best partner, with k)
  protein_ids.txt      all query and partner sequence IDs ('<species>__<gene>'), for extracting proteins
  cds_ids/<species>.txt  gene IDs per species, for extracting coding sequences
"""
import argparse
import collections
import csv
import os

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--scan", nargs="+", required=True, help="DIAMOND table(s); give the first scan, then any masking-off rescan")
ap.add_argument("--accessions", required=True, help="data/txptome_accessions.csv (Species, Family columns)")
ap.add_argument("--family-override", nargs="*", default=["Phrurolithus_coreanus=Phrurolithidae"],
                help="SPECIES=Family corrections to the accession table")
ap.add_argument("--min-pident", type=float, default=98.0)
ap.add_argument("--min-aligned-aa", type=int, default=100)
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()

fam = {r["Species"]: r["Family"] for r in csv.DictReader(open(args.accessions, encoding="utf-8-sig"))}
for kv in args.family_override:
    sp, f = kv.split("=")
    fam[sp] = f

cols = ["q", "s", "pid", "len", "qcov", "scov", "e", "bits"]
df = pd.concat([pd.read_csv(p, sep="\t", header=None, names=cols) for p in args.scan])
df["qsp"] = df.q.str.split("__").str[0]
df["ssp"] = df.s.str.split("__").str[0]
df = df[df.qsp != df.ssp]
df["qfam"] = df.qsp.map(fam)
df["sfam"] = df.ssp.map(fam)
x = df[(df.qfam != df.sfam) & (df.pid >= args.min_pident) & (df.len >= args.min_aligned_aa)]
g = x.groupby("q").agg(nsp=("ssp", "nunique"))
g["k"] = g.nsp + 1
best = (x.sort_values(["pid", "len"], ascending=False).drop_duplicates("q")
        .merge(g[["k"]], left_on="q", right_index=True))
print(f"{len(best)} suspect genes with a cross-family partner (>= {args.min_pident}% identity, >= {args.min_aligned_aa} aa)")

os.makedirs(os.path.join(args.out_dir, "cds_ids"), exist_ok=True)
best.to_csv(os.path.join(args.out_dir, "candidates.tsv"), sep="\t", index=False)
ids = set(best.q) | set(best.s)
open(os.path.join(args.out_dir, "protein_ids.txt"), "w").write("\n".join(sorted(ids)) + "\n")
by = collections.defaultdict(list)
for i in ids:
    sp, gid = i.split("__")
    by[sp].append(gid)
for sp, l in by.items():
    open(os.path.join(args.out_dir, "cds_ids", f"{sp}.txt"), "w").write("\n".join(l) + "\n")
print(f"{len(ids)} sequences to extract from {len(by)} species")
