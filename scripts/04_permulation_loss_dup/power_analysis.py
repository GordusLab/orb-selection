"""Power of the odds ratio permulation test, using synthetic genes.

Question: if genes really tracked orb-weaving, how often would the test (permulation p-values plus
Benjamini-Hochberg correction) flag them? For each scenario, synthetic genes are generated whose
loss (or duplication) pattern follows the true foreground labeling to a chosen degree, tested
against every permulation exactly like a real HOG, and then mixed into the real HOGs' p-values to
see whether BH would call them when 1, 2, 5 or 20 such genes exist.

Scenarios are "effect vs other" pairs of event probabilities:
  * loss: the gene is absent with probability `effect` in non-orb-weavers (background) and `other`
    in orb-weavers (foreground), so loss is more likely in the background;
  * duplication: the gene is duplicated with probability `effect` in orb-weavers and `other` in
    non-orb-weavers.
"1.0 vs 0.0" is a gene that follows the labeling perfectly; "0.5 vs 0.5" is a negative control.

Two corrections are reported, because the two directions of a test are two tails of one comparison:
  * stacked: both tails of each test corrected together (the two-sided analysis);
  * directional: only the pre-specified tail (loss in non-orb-weavers, duplication in orb-weavers),
    corrected within that tail.

Caveat: the synthetic genes have independent species-by-species noise, whereas real genes tend to
follow clades, which a clade-structured permulation can mimic more easily. The detection rates are
therefore best read as an upper bound on real power.

Usage (from anywhere in the repo):
    python scripts/04_permulation_loss_dup/power_analysis.py \
        --permulations-csv data/perms_tip_values_1M.csv --out-csv results/power_analysis_1M.csv
"""

import argparse
import os
import sys
import time

import matplotlib

matplotlib.use("Agg")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_HERE, os.path.join(_HERE, os.pardir, os.pardir, "src")]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import odds_ratio_test as ort  # noqa: E402
from odds_ratio_test_helpers import (  # noqa: E402
    _resolve_repo_path,
    benjamini_hochberg,
    calculate_odds,
    calculate_odds_batch,
    load_permulation_matrices,
)

SCENARIOS = ((1.0, 0.0), (0.9, 0.1), (0.8, 0.2), (0.7, 0.3), (0.6, 0.4), (0.5, 0.5))
SPIKE_COUNTS = (1, 2, 5, 20)
# test name -> (tail where the effect shows up, group in which the event is more frequent)
TESTS = {"loss": ("bg", "background"), "dup": ("fg", "foreground")}


def tail_p_values(mat, busco, observed, foreground, background, tol, chunk_size=50):
    """Add-one permulation p-values (ties inclusive) in both directions for each row of `mat`.

    Returns (p_fg, p_bg): "fg" = log odds ratio at least as large as the permulated values (event
    more likely in the foreground), "bg" = at least as small (more likely in the background).
    """
    n_perm = foreground.shape[1]
    p_fg = np.empty(mat.shape[0])
    p_bg = np.empty(mat.shape[0])
    for start in range(0, mat.shape[0], chunk_size):
        stop = min(start + chunk_size, mat.shape[0])
        null = calculate_odds_batch(foreground, background, mat[start:stop], busco)
        obs = observed[start:stop, None]
        p_fg[start:stop] = ((null >= obs - tol).sum(axis=1) + 1) / (n_perm + 1)
        p_bg[start:stop] = ((null <= obs + tol).sum(axis=1) + 1) / (n_perm + 1)
    return p_fg, p_bg


def fraction_called(real_p, spiked_p, other_tail_p, q_level=0.05):
    """Fraction of the spiked genes that BH calls significant when added to the real p-values."""
    stacked = np.concatenate([real_p, spiked_p, other_tail_p])
    q = benjamini_hochberg(stacked)
    return float(np.mean(q[len(real_p) : len(real_p) + len(spiked_p)] < q_level))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--permulations-csv", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--n-genes", type=int, default=300, help="Synthetic genes per scenario")
    parser.add_argument("--reps", type=int, default=100, help="Random draws per spike count")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--min-occ", type=int, default=30)
    parser.add_argument("--max-occ", type=int, default=88, help="Upper occupancy bound for loss")
    args = parser.parse_args()

    root = _resolve_repo_path(".")
    true_odds = ort.OddsRatioResults(
        genecount_csv=os.path.join(root, "data/N5.GeneCount.tsv"),
        hog_node_genes_tsv=os.path.join(root, "data/N5.tsv"),
        time="power analysis",
        foreground_list_filename=os.path.join(root, "data/orbweavers-list.txt"),
        background_list_filename=None,
        buscos_filename=os.path.join(root, "data/buscos.csv"),
    )
    included = (true_odds.foreground_bool_arr + true_odds.background_bool_arr).astype(bool).nonzero()[0]
    foreground, background = load_permulation_matrices(
        _resolve_repo_path(args.permulations_csv), true_odds.all_species_arr, included
    )
    n_perm = foreground.shape[1]
    floor = 1 / (n_perm + 1)
    tol = ort.TIE_TOLERANCE
    is_foreground = true_odds.foreground_bool_arr.ravel().astype(bool)
    occ = true_odds.total_occ_arr.ravel()
    rng = np.random.default_rng(args.seed)
    rows_out = []
    t0 = time.time()

    for test, (effect_tail, event_group) in TESTS.items():
        mat_all = getattr(true_odds, f"{test}_bool_mat")
        busco = getattr(true_odds, f"{test}_busco_arr")
        observed_all = getattr(true_odds, f"{test}_lor_arr").ravel()
        in_range = occ >= args.min_occ
        if test == "loss":
            in_range &= occ <= args.max_occ
        rows = np.where(in_range)[0]

        # the real HOGs' p-values under this permulation set
        real_fg, real_bg = tail_p_values(mat_all[rows], busco, observed_all[rows], foreground, background, tol)
        real_stacked = np.concatenate([real_fg, real_bg])
        real_effect_tail = real_bg if effect_tail == "bg" else real_fg
        print(
            f"[{test}] {len(rows)} real HOGs, {n_perm} permulations ({time.time() - t0:.0f}s); smallest p "
            f"{real_stacked.min():.2e}; BH needs {int(np.ceil(floor * len(real_stacked) / 0.05))} genes at the "
            f"floor ({floor:.1e}) with both tails stacked, {int(np.ceil(floor * len(rows) / 0.05))} with one tail",
            flush=True,
        )

        in_event_group = is_foreground if event_group == "foreground" else ~is_foreground
        for effect, other in SCENARIOS:
            prob = np.where(in_event_group, effect, other)
            events = (rng.random((args.n_genes, len(is_foreground))) < prob).astype(float)
            observed = np.array(
                [
                    calculate_odds(
                        true_odds.foreground_bool_arr, true_odds.background_bool_arr, events[i : i + 1], busco_arr=busco
                    ).ravel()[0]
                    for i in range(args.n_genes)
                ]
            )
            p_fg, p_bg = tail_p_values(events, busco, observed, foreground, background, tol)
            p_effect = p_bg if effect_tail == "bg" else p_fg

            row = {
                "test": test,
                "effect": f"{effect:.1f} vs {other:.1f}",
                "median_p": float(np.median(p_effect)),
                "frac_p_le_1e-4": float(np.mean(p_effect <= 1e-4)),
                "frac_p_le_1e-5": float(np.mean(p_effect <= 1e-5)),
                "frac_at_floor": float(np.mean(p_effect <= floor * 1.001)),
            }
            for k in SPIKE_COUNTS:
                stacked_hits, directional_hits = [], []
                for _ in range(args.reps):
                    idx = rng.choice(args.n_genes, k, replace=False)
                    stacked_hits.append(fraction_called(real_stacked, p_effect[idx], np.ones(k)))
                    directional_hits.append(fraction_called(real_effect_tail, p_effect[idx], np.empty(0)))
                row[f"BH_stacked_k{k}"] = float(np.mean(stacked_hits))
                row[f"BH_directional_k{k}"] = float(np.mean(directional_hits))
            rows_out.append(row)
            print(f"  {test} {row['effect']}: done ({time.time() - t0:.0f}s)", flush=True)

    result = pd.DataFrame(rows_out)
    out_path = _resolve_repo_path(args.out_csv)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    result.to_csv(out_path, index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(f"\nWrote {out_path}\n")
        print(result.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
