"""Run the odds ratio permulation test from the command line (e.g. as a SLURM job).

Wraps odds_ratio_test.odds_ratio_test() with the same defaults as the notebook. Paths that are not
absolute are resolved from the repository root.

Example (from anywhere inside the repo):
    python scripts/04_permulation_loss_dup/run_odds_ratio_test.py \
        --permulations-csv data/perms_tip_values_100k.csv --dir-suffix 100k
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")  # no display on cluster nodes

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_HERE, os.path.join(_HERE, os.pardir, os.pardir, "src")]

import odds_ratio_test as ort  # noqa: E402
from odds_ratio_test_helpers import _resolve_repo_path  # noqa: E402


def count_csv_rows(path):
    """Number of data rows (excluding the header) in a CSV."""
    with open(path, encoding="utf-8") as file_obj:
        return sum(1 for _ in file_obj) - 1


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--permulations-csv", default="data/perms_tip_values.csv")
    parser.add_argument(
        "--reps",
        type=int,
        default=None,
        help="Number of permulations to use (default: every row of the CSV)",
    )
    parser.add_argument("--min-occ", type=int, default=30)
    parser.add_argument("--max-occ", type=int, default=88, help="Upper occupancy bound for loss")
    parser.add_argument("--alpha", type=float, default=0.05, help="Two-sided level")
    parser.add_argument("--results-dir", default="results/odds_ratio_test")
    parser.add_argument("--dir-suffix", default=None)
    parser.add_argument("--species-of-interest", default="Uloborus_diversus")
    parser.add_argument("--fg-name", default="orbweavers")
    parser.add_argument("--bg-name", default="non-orbweavers")
    parser.add_argument(
        "--no-pickle",
        action="store_true",
        help="Do not save results.pkl (it includes every permulated tip assignment)",
    )
    args = parser.parse_args()

    n_rows = count_csv_rows(_resolve_repo_path(args.permulations_csv))
    reps = args.reps if args.reps is not None else n_rows

    ort.odds_ratio_test(
        min_occ=args.min_occ,
        max_occ=args.max_occ,
        alpha=args.alpha,
        permulation_reps=reps,
        permulations_tip_values_csv=args.permulations_csv,
        species_of_interest=args.species_of_interest,
        results_dir=args.results_dir,
        fg_name=args.fg_name,
        bg_name=args.bg_name,
        save_pickle=not args.no_pickle,
        dir_suffix=args.dir_suffix,
    )
    print("DONE")


if __name__ == "__main__":
    main()
