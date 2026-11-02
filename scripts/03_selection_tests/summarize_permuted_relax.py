#!/usr/bin/env python3
"""
Compare real RELAX results to the permuted-labeling null from run_relax_permuted.sh.

For each permuted labeling, computed on the sampled HOGs it completed:
  - fraction of HOGs with RELAX p <= alpha
  - fraction of HOGs with k > 1 (all HOGs, significant or not)
  - median log2(k)
The same statistics are computed for the real labeling on the same HOGs. The
empirical p-value is (1 + #labelings with stat >= real) / (1 + #labelings), so
its floor is 1 / (1 + n_labelings). Only labelings that finished every sampled
HOG are used, so each statistic is over the same set of genes.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hyphy_results_parser import HyphyResultsManager, HyphyResult


def stats_for(df, alpha):
    return {
        "frac_sig": (df["p_value"].astype(float) <= alpha).mean(),
        "frac_k_gt_1": (df["k"].astype(float) > 1).mean(),
        "median_log2_k": np.log2(df["k"].astype(float).clip(1e-3, 1e3)).median(),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perm-root", required=True, help="OUT_ROOT of run_relax_permuted.sh")
    ap.add_argument("--hogs", required=True, help="hogs.txt from make_permuted_relax_manifest.py")
    ap.add_argument("--real", default="../../results/hyphy_results_cache/relax_results.pkl")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", default=None, help="optional CSV of per-labeling statistics")
    args = ap.parse_args()

    hogs = [line.strip() for line in open(args.hogs) if line.strip()]
    real = HyphyResult.load_from_pickle(args.real).results_df.loc[hogs]
    real_stats = stats_for(real, args.alpha)

    rows = []
    for perm_id in sorted(os.listdir(args.perm_root)):
        d = os.path.join(args.perm_root, perm_id)
        if not os.path.isdir(d):
            continue
        res = HyphyResultsManager().load_relax_from_json(d, name=perm_id).results_df
        if not set(hogs) <= set(res.index):
            continue
        rows.append({"perm_id": perm_id, **stats_for(res.loc[hogs], args.alpha)})

    null = pd.DataFrame(rows)
    print(f"{len(hogs)} HOGs; {len(null)} complete labelings\n")
    for col, real_val in real_stats.items():
        n_ge = int((null[col] >= real_val).sum())
        n_le = int((null[col] <= real_val).sum())
        print(f"{col:14s} real={real_val:7.3f}  null mean={null[col].mean():7.3f} "
              f"[{null[col].min():.3f}, {null[col].max():.3f}]  "
              f"p(null>=real)={(1 + n_ge) / (1 + len(null)):.3f}  p(null<=real)={(1 + n_le) / (1 + len(null)):.3f}")
    if args.out:
        null.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
