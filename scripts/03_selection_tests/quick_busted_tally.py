#!/usr/bin/env python3
"""Quick tally of BUSTED-PH hits from completed JSONs, without pandas/numpy.

A "hit" (per hyphy_results_parser.py's load_busted_ph_from_json): test p <= 0.05,
background p > 0.05, shared p <= 0.05.

Usage: python3 quick_busted_tally.py <glob-pattern>
Example: python3 quick_busted_tally.py \
    "/scratch4/agordus1/crunnel2/hyphy_wd_260929/*/*_BUSTED-PH_orb_fg.json"
"""
import glob
import json
import sys


def first_present(mapping, *keys):
    for key in keys:
        if key in mapping:
            return mapping[key]
    raise KeyError(f"None of the expected keys were found: {', '.join(keys)}")


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <glob-pattern>")

    files = sorted(glob.glob(sys.argv[1]))
    n_complete = 0
    n_hit = 0
    n_error = 0

    for path in files:
        try:
            with open(path) as f:
                data = json.load(f)
            test_p = float(data["test results"]["p-value"])
            bg_p = float(first_present(
                data, "test results background", "Background selection test results",
            )["p-value"])
            shared_p = float(first_present(
                data, "test results shared distributions", "Comparative selection test results",
            )["p-value"])
        except Exception as e:
            n_error += 1
            print(f"SKIP {path}: {e}", file=sys.stderr)
            continue

        n_complete += 1
        if test_p <= 0.05 and bg_p > 0.05 and shared_p <= 0.05:
            n_hit += 1
            print(f"HIT\t{path}\t{test_p}\t{bg_p}\t{shared_p}")

    print(
        f"\n{n_complete} completed JSONs matched, {n_error} skipped, "
        f"{n_hit} hits ({100 * n_hit / n_complete:.1f}%)" if n_complete else
        f"\n0 completed JSONs matched, {n_error} skipped",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
