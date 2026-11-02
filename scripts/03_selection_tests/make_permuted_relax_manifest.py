#!/usr/bin/env python3
"""
Build the task manifest for run_relax_permuted.sh.

Null-calibration design: a fixed random subset of HOGs (--n-hogs) is run under
a set of permuted phenotype labelings (--n-perms), so every labeling is applied
to the SAME HOGs. This keeps the labeling-level correlation between genes (a
labeling that puts the foreground in a fast clade shifts every gene), which is
what the genome-wide summary statistics (fraction significant, fraction k > 1)
need to be compared against. One random labeling per HOG would average that
variance away.

Each permuted row of --perms is 1 for tips assigned "orb-weaver" (same count as
the real phenotype), so the RELAX foreground (non-orb-weavers in the real
analysis) is the tips with value 0.

Writes <outdir>/manifest.tsv with columns: HOG, perm_id, regex. The regex has
the same form prep_for_hyphy2.sh passes to label-tree.bf.
"""

import argparse
import json
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perms", default="data/perms_tip_values_100k.csv", help="perm_id + one 0/1 column per tip")
    ap.add_argument("--hog-lists", nargs="+", default=["data/hyphy_tiers/tier_A.txt"],
                    help="HOG lists to sample from (default: runtime tier A only, the cheapest)")
    ap.add_argument("--n-hogs", type=int, default=60)
    ap.add_argument("--n-perms", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--runtime-json", default=None,
                    help="path template of each HOG's real RELAX JSON, e.g. WD/{HOG}/{HOG}_RELAX.json; "
                         "if given, tasks are ordered longest real runtime first so the array's tail is short")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    hogs = sorted({line.strip() for f in args.hog_lists for line in open(f) if line.strip()})
    hogs = sorted(rng.choice(hogs, size=min(args.n_hogs, len(hogs)), replace=False))

    perms = pd.read_csv(args.perms)
    perms = perms.iloc[rng.choice(len(perms), size=min(args.n_perms, len(perms)), replace=False)]
    tips = list(perms.columns[1:])

    regex_by_perm = {}
    for _, row in perms.iterrows():
        fg_tips = [t for t in tips if row[t] == 0]
        regex_by_perm[row["perm_id"]] = "^(" + "|".join(fg_tips) + ")\\|"

    os.makedirs(args.outdir, exist_ok=True)
    rows = [(h, p, regex_by_perm[p]) for p in regex_by_perm for h in hogs]
    if args.runtime_json:
        hours = {h: json.load(open(args.runtime_json.format(HOG=h)))["timers"]["Overall"]["timer"] / 3600
                 for h in hogs}
        rows.sort(key=lambda r: -hours[r[0]])
        print(f"real RELAX runtime of sampled HOGs (h): mean {np.mean(list(hours.values())):.1f}, "
              f"max {max(hours.values()):.1f}; tasks ordered longest first")
    pd.DataFrame(rows, columns=["HOG", "perm_id", "regex"]).to_csv(
        os.path.join(args.outdir, "manifest.tsv"), sep="\t", index=False, header=False
    )
    pd.Series(hogs).to_csv(os.path.join(args.outdir, "hogs.txt"), index=False, header=False)
    print(f"{len(hogs)} HOGs x {len(regex_by_perm)} labelings = {len(rows)} tasks "
          f"-> {args.outdir}/manifest.tsv (use --array=1-{len(rows)})")


if __name__ == "__main__":
    main()
