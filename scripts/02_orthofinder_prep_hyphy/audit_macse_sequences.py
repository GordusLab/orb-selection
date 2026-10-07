#!/usr/bin/env python3
"""Audit MACSE NT alignments for sequences that might be filtered after alignment.

For every <WD>/<HOG>/macse/<HOG>_NT.fasta, reports per sequence:
  - ungapped length and its fraction of the gene's median ungapped length
  - number of frameshift marks ("!")
  - number of internal in-frame stop codons (TAA/TAG/TGA, excluding the final codon)
and flags sequences below --min-frac of the median, or with frameshifts / internal stops.
Flagged sequences are cross-referenced with the orbweaver list and Uloborus_diversus so
you can see which genes would lose the taxa that make them testable.

Usage:
  audit_macse_sequences.py <WD> <hog_list.txt> <orbweavers-list.txt> <out_prefix> [--min-frac 0.5] [--max-internal-stops 0] [--udiv Uloborus_diversus]

Writes <out_prefix>.sequences.tsv, <out_prefix>.genes.tsv and prints a summary.
Nothing is modified or removed.
"""
import argparse
import csv
import os
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor

STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    names, seqs = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                names.append(line[1:])
                seqs.append([])
            else:
                seqs[-1].append(line)
    return names, ["".join(s) for s in seqs]


def internal_stops(seq):
    """In-frame stop codons, excluding the last codon that has any sequence."""
    seq = seq.upper()
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    last = -1
    for i, c in enumerate(codons):
        if c != "---":
            last = i
    return sum(1 for i, c in enumerate(codons) if i < last and c in STOPS)


def audit_gene(args):
    wd, hog = args
    path = os.path.join(wd, hog, "macse", f"{hog}_NT.fasta")
    if not os.path.isfile(path):
        return hog, None
    names, seqs = read_fasta(path)
    if not seqs:
        return hog, []
    lens = [sum(1 for ch in s if ch not in "-!") for s in seqs]
    rows = []
    for name, s, n in zip(names, seqs, lens):
        rows.append(
            {
                "hog": hog,
                "sequence": name,
                "species": name.split("|")[0],
                "aligned_len": len(s),
                "ungapped_len": n,
                "n_frameshift": s.count("!"),
                "n_internal_stops": internal_stops(s),
            }
        )
    med = statistics.median(lens)
    for r in rows:
        r["frac_of_median"] = round(r["ungapped_len"] / med, 4) if med else 0.0
    return hog, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wd")
    ap.add_argument("hog_list")
    ap.add_argument("orb_list")
    ap.add_argument("out_prefix")
    ap.add_argument("--min-frac", type=float, default=0.5, help="flag sequences below this fraction of the median ungapped length")
    ap.add_argument("--max-internal-stops", type=int, default=0, help="flag sequences with more internal stops than this")
    ap.add_argument("--udiv", default="Uloborus_diversus")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    with open(a.hog_list) as fh:
        hogs = [l.strip() for l in fh if l.strip()]
    with open(a.orb_list) as fh:
        orbs = {l.strip().replace("\r", "") for l in fh if l.strip()}

    seq_rows, gene_rows, missing = [], [], []
    with ProcessPoolExecutor(max_workers=a.threads) as ex:
        results = ex.map(audit_gene, [(a.wd, h) for h in hogs], chunksize=8)
        for hog, rows in results:
            if rows is None:
                missing.append(hog)
                continue
            if not rows:
                continue
            for r in rows:
                r["is_orb"] = int(r["species"] in orbs)
                r["is_udiv"] = int(r["species"] == a.udiv)
                r["short"] = int(r["frac_of_median"] < a.min_frac)
                r["has_fs"] = int(r["n_frameshift"] > 0)
                r["has_stops"] = int(r["n_internal_stops"] > a.max_internal_stops)
                r["flagged"] = int(r["short"] or r["has_fs"] or r["has_stops"])
            seq_rows.extend(rows)

            kept = [r for r in rows if not r["flagged"]]
            gene_rows.append(
                {
                    "hog": hog,
                    "n_seqs": len(rows),
                    "median_ungapped_len": statistics.median(r["ungapped_len"] for r in rows),
                    "n_short": sum(r["short"] for r in rows),
                    "n_frameshift_seqs": sum(r["has_fs"] for r in rows),
                    "n_stop_seqs": sum(r["has_stops"] for r in rows),
                    "n_flagged": sum(r["flagged"] for r in rows),
                    "orb_total": sum(r["is_orb"] for r in rows),
                    "orb_flagged": sum(r["is_orb"] and r["flagged"] for r in rows),
                    "orb_remaining": sum(r["is_orb"] for r in kept),
                    "nonorb_remaining": sum(not r["is_orb"] for r in kept),
                    "udiv_flagged": int(any(r["is_udiv"] and r["flagged"] for r in rows)),
                    "udiv_present": int(any(r["is_udiv"] for r in rows)),
                }
            )

    seq_fields = [
        "hog", "sequence", "species", "is_orb", "is_udiv", "aligned_len", "ungapped_len",
        "frac_of_median", "n_frameshift", "n_internal_stops", "short", "has_fs", "has_stops", "flagged",
    ]
    with open(f"{a.out_prefix}.sequences.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=seq_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(seq_rows)
    gene_fields = list(gene_rows[0].keys()) if gene_rows else []
    with open(f"{a.out_prefix}.genes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=gene_fields, delimiter="\t")
        w.writeheader()
        w.writerows(gene_rows)

    n_g = len(gene_rows)
    n_s = len(seq_rows)
    flagged = [r for r in seq_rows if r["flagged"]]
    print(f"genes audited: {n_g}  (missing macse NT file: {len(missing)})")
    print(f"sequences: {n_s}  flagged: {len(flagged)} ({100 * len(flagged) / max(n_s, 1):.2f}%)")
    print(f"  short (<{a.min_frac} of median): {sum(r['short'] for r in seq_rows)}")
    print(f"  with frameshift '!': {sum(r['has_fs'] for r in seq_rows)}")
    print(f"  with >{a.max_internal_stops} internal stops: {sum(r['has_stops'] for r in seq_rows)}")
    print(f"genes with >=1 flagged sequence: {sum(g['n_flagged'] > 0 for g in gene_rows)}")
    print(f"genes where U.div sequence is flagged: {sum(g['udiv_flagged'] for g in gene_rows)}")
    print(f"genes that would have 0 orbweavers left: {sum(g['orb_total'] > 0 and g['orb_remaining'] == 0 for g in gene_rows)}")
    print(f"genes that would have 0 non-orbweavers left: {sum(g['nonorb_remaining'] == 0 for g in gene_rows)}")
    if missing:
        print(f"first missing: {missing[:5]}", file=sys.stderr)


if __name__ == "__main__":
    main()
