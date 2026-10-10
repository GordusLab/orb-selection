#!/usr/bin/env python3
"""Quick tally of RELAX results from completed JSONs, without pandas/numpy.

Classification mirrors hyphy_results_parser.py's load_relax_from_json:
p <= 0.05 & k < 1  -> relaxed
p <= 0.05 & k > 1  -> intensified
otherwise          -> not significant

Usage: python3 quick_relax_tally.py <glob-pattern>
Example: python3 quick_relax_tally.py \
    "/scratch4/agordus1/crunnel2/hyphy_wd_260929/*/*_RELAX.json"
"""
import glob
import json
import sys
from collections import Counter


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <glob-pattern>")

    files = sorted(glob.glob(sys.argv[1]))
    n_error = 0
    counts = Counter()

    for path in files:
        try:
            with open(path) as f:
                data = json.load(f)
            pval = float(data["test results"]["p-value"])
            k = float(data["test results"]["relaxation or intensification parameter"])
        except Exception as e:
            n_error += 1
            print(f"SKIP {path}: {e}", file=sys.stderr)
            continue

        if pval <= 0.05 and k < 1:
            result = "relaxed"
        elif pval <= 0.05 and k > 1:
            result = "intensified"
        else:
            result = "not significant"

        counts[result] += 1
        print(f"{result}\t{path}\t{pval}\t{k}")

    n_complete = sum(counts.values())
    print(
        f"\n{n_complete} completed JSONs matched, {n_error} skipped: "
        + ", ".join(f"{k}={v}" for k, v in counts.most_common()),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
