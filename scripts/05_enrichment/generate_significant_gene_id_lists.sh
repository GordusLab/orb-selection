#!/usr/bin/env bash
set -euo pipefail

# Generate significant_gene_id_lists outputs from pickle files.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
SCRIPT="$REPO_ROOT/scripts/05_enrichment/get_gene_id_lists.py"
OUT_BASE="$REPO_ROOT/results/significant_gene_id_lists"

LOSS_PKL="$REPO_ROOT/results/odds_ratio_test/Results_Mar25/Run1_Loss_LT_50-95_10000x/loss_occ50-95_less.pkl"
DUP_PKL="$REPO_ROOT/results/odds_ratio_test/Results_Mar25/Run2_Dup_RT_50-max_10000x/duplication_occ50-98_greater.pkl"
RELAX_PKL="$REPO_ROOT/results/hyphy_results_cache/relax_results.pkl"
BUSTED_PKL="$REPO_ROOT/results/hyphy_results_cache/busted_ph_results.pkl"
BUSTED_REV_PKL="$REPO_ROOT/results/hyphy_results_cache/busted_ph_rev_results.pkl"

mkdir -p "$OUT_BASE/loss" "$OUT_BASE/duplication" "$OUT_BASE/relax" "$OUT_BASE/busted_ph" "$OUT_BASE/busted_ph_rev"

echo "Generating odds-ratio/permulation ID lists..."
"$PYTHON_BIN" "$SCRIPT" "$LOSS_PKL" --tail left  --hits-file "$OUT_BASE/loss/loss_nonorb.txt" --universe-file "$OUT_BASE/loss/loss_universe.txt"
"$PYTHON_BIN" "$SCRIPT" "$LOSS_PKL" --tail right --hits-file "$OUT_BASE/loss/loss_orb.txt"

"$PYTHON_BIN" "$SCRIPT" "$DUP_PKL" --tail left  --hits-file "$OUT_BASE/duplication/duplication_nonorb.txt" --universe-file "$OUT_BASE/duplication/duplication_universe.txt"
"$PYTHON_BIN" "$SCRIPT" "$DUP_PKL" --tail right --hits-file "$OUT_BASE/duplication/duplication_orb.txt"

echo "Generating HyPhy ID lists..."
"$PYTHON_BIN" "$SCRIPT" "$RELAX_PKL" --hits-file "$OUT_BASE/relax/relax_all_hits.txt" --universe-file "$OUT_BASE/relax/relax_universe.txt"
"$PYTHON_BIN" "$SCRIPT" "$RELAX_PKL" --relax-result relaxed --hits-file "$OUT_BASE/relax/relax_relaxed_hits.txt"
"$PYTHON_BIN" "$SCRIPT" "$RELAX_PKL" --relax-result intensified --hits-file "$OUT_BASE/relax/relax_intensified_hits.txt"

"$PYTHON_BIN" "$SCRIPT" "$BUSTED_PKL" --hits-file "$OUT_BASE/busted_ph/busted_ph_hits.txt" --universe-file "$OUT_BASE/busted_ph/busted_ph_universe.txt"
"$PYTHON_BIN" "$SCRIPT" "$BUSTED_REV_PKL" --hits-file "$OUT_BASE/busted_ph_rev/busted_ph_rev_hits.txt" --universe-file "$OUT_BASE/busted_ph_rev/busted_ph_rev_universe.txt"

echo "Done. Regenerated significant_gene_id_lists files under:"
echo "  $OUT_BASE"
