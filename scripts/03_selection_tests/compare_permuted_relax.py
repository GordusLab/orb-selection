#!/usr/bin/env python3
"""
Compare the real RELAX results with the permuted-label null, from a table of per-run results.

Two steps (see run_relax_permuted.sh and make_permuted_relax_manifest.py for the runs):

  extract  read <PERM_ROOT>/<perm_id>/<HOG>_RELAX.json (finished runs only) into a TSV of
           perm_id, HOG, p_value, k. Uses only the standard library, so it can run on the cluster.
  compare  for each permuted labeling, compute statistics over the sampled HOGs and compare
           them with the real labeling on the same HOGs. Statistics: fraction of HOGs with p <= alpha,
           fraction with k > 1, median log2 k, and the fractions that are significant and intensified
           (k > 1) or relaxed (k < 1). Empirical p = (1 + #labelings with stat >= real) / (1 + #labelings).
           Only labelings with every sampled HOG finished are used.

Usage:
  compare_permuted_relax.py extract --perm-root DIR --out perm_results.tsv
  compare_permuted_relax.py compare --results perm_results.tsv --hogs hogs.txt [--real PKL] [--alpha 0.05] [--out stats.csv]
"""

import argparse
import csv
import glob
import json
import os
import sys


def extract(perm_root, out):
    rows = []
    for path in glob.glob(os.path.join(perm_root, "*", "*_RELAX.json")):
        try:
            test = json.load(open(path)).get("test results")
        except (OSError, ValueError):
            continue
        if not test:
            continue
        rows.append((os.path.basename(os.path.dirname(path)), os.path.basename(path)[: -len("_RELAX.json")],
                     test["p-value"], test["relaxation or intensification parameter"]))
    with open(out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["perm_id", "HOG", "p_value", "k"])
        writer.writerows(sorted(rows))
    print(f"{len(rows)} finished runs written to {out}")


def stats_for(df, alpha):
    import numpy as np

    k = df["k"].astype(float)
    p = df["p_value"].astype(float)
    return {
        "frac_sig": (p <= alpha).mean(),
        "frac_k_gt_1": (k > 1).mean(),
        "median_log2_k": float(np.log2(k.clip(1e-3, 1e3)).median()),
        "frac_sig_intensified": ((p <= alpha) & (k > 1)).mean(),
        "frac_sig_relaxed": ((p <= alpha) & (k < 1)).mean(),
    }


def compare(results, hogs_file, real_pkl, alpha, out):
    import pandas as pd

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hyphy_results_parser import HyphyResult

    hogs = [line.strip() for line in open(hogs_file) if line.strip()]
    runs = pd.read_csv(results, sep="\t")
    runs = runs[runs["HOG"].isin(hogs)]
    complete = runs.groupby("perm_id")["HOG"].nunique()
    complete = complete[complete == len(hogs)].index
    runs = runs[runs["perm_id"].isin(complete)]
    real = HyphyResult.load_from_pickle(real_pkl).results_df.loc[hogs]
    real_stats = stats_for(real, alpha)
    null = runs.groupby("perm_id").apply(lambda d: pd.Series(stats_for(d, alpha)))
    print(f"{len(hogs)} HOGs; {len(null)} complete labelings\n")
    print(f"{'statistic':22s}{'real':>8s}{'null mean':>11s}{'2.5%':>8s}{'97.5%':>8s}  p(null >= real)")
    for col, value in real_stats.items():
        n_ge = int((null[col] >= value).sum())
        print(f"{col:22s}{value:8.3f}{null[col].mean():11.3f}{null[col].quantile(.025):8.3f}"
              f"{null[col].quantile(.975):8.3f}  {(1 + n_ge) / (1 + len(null)):.3f}")
    if out:
        null.to_csv(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--perm-root", required=True)
    e.add_argument("--out", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--results", required=True)
    c.add_argument("--hogs", required=True)
    c.add_argument("--real", default="../../results/hyphy_results_cache/relax_results.pkl")
    c.add_argument("--alpha", type=float, default=0.05)
    c.add_argument("--out", default=None, help="optional CSV of per-labeling statistics")
    args = ap.parse_args()
    if args.command == "extract":
        extract(args.perm_root, args.out)
    else:
        compare(args.results, args.hogs, args.real, args.alpha, args.out)


if __name__ == "__main__":
    main()
