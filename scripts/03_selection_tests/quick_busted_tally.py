#!/usr/bin/env python3
"""Quick tally of BUSTED-PH hits from completed JSONs, without pandas/numpy.

A "hit" is read directly from HyPhy's own "Summary" sentence rather than
re-deriving it from the three p-values (which has version-dependent key
names across HyPhy builds, e.g. "test results background" vs "Background
selection test results"): a hit is a Summary containing "associated with
the trait" without a "**no**" right before it.

Usage: python3 quick_busted_tally.py <glob-pattern>
Example: python3 quick_busted_tally.py \
    "/scratch4/agordus1/crunnel2/hyphy_wd_260929/*/*_BUSTED-PH_orb_fg.json"
"""
import glob
import json
import sys


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
            summary = data["BUSTED-PH"]["Summary"]
        except Exception:
            # Most skips are genes still mid-run (incomplete JSON); not worth printing per-file.
            n_error += 1
            continue

        n_complete += 1
        if "associated with the trait" in summary and "**no**" not in summary:
            n_hit += 1
            print(f"HIT\t{path}")

    print(
        f"\n{n_complete} completed JSONs matched, {n_error} skipped, "
        f"{n_hit} hits ({100 * n_hit / n_complete:.1f}%)" if n_complete else
        f"\n0 completed JSONs matched, {n_error} skipped",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
