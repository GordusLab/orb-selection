"""Combine permulation batch CSVs into one tip-values file and sanity-check the result.

Usage (from the repo root):
    python scripts/04_permulation_loss_dup/combine_permulation_batches.py \
        --include data/perms_tip_values.csv \
        --out data/perms_tip_values_100k.csv
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches-dir", default="data/perms_batches")
    parser.add_argument(
        "--include",
        nargs="*",
        default=[],
        help="Existing tip-value CSVs to include (e.g. the original 10,000 permulations)",
    )
    parser.add_argument("--out", default="data/perms_tip_values_100k.csv")
    parser.add_argument(
        "--max-batch",
        type=int,
        default=None,
        help="Only use batches with an id up to this number (batches are named batch_<id>_tip_values.csv)",
    )
    args = parser.parse_args()

    def batch_id(path):
        return int(re.search(r"batch_(\d+)_tip_values", os.path.basename(path)).group(1))

    batch_paths = sorted(
        glob.glob(os.path.join(args.batches_dir, "batch_*_tip_values.csv")), key=batch_id
    )
    if args.max_batch is not None:
        batch_paths = [p for p in batch_paths if batch_id(p) <= args.max_batch]
    paths = list(args.include) + batch_paths
    if not paths:
        sys.exit("No input CSVs found.")

    frames = []
    for path in paths:
        df = pd.read_csv(path)
        print(f"{path}: {len(df)} permulations")
        frames.append(df.drop(columns="perm_id"))

    species = list(frames[0].columns)
    for path, df in zip(paths, frames):
        if set(df.columns) != set(species):
            sys.exit(f"{path} has different species columns than {paths[0]}")

    combined = pd.concat([df[species] for df in frames], ignore_index=True)

    fg_counts = combined.sum(axis=1)
    n_unique = len(combined.drop_duplicates())
    print(
        f"\nTotal permulations: {len(combined)}\n"
        f"Species: {len(species)}\n"
        f"Foreground count per permulation: min {fg_counts.min():g}, max {fg_counts.max():g}\n"
        f"Unique tip-value maps: {n_unique} ({len(combined) - n_unique} repeats)"
    )
    if fg_counts.min() != fg_counts.max():
        print("WARNING: foreground count varies between permulations.")

    combined.insert(0, "perm_id", range(1, len(combined) + 1))
    combined.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
