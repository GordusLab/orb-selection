#!/usr/bin/env bash
# Run several permulation batches as parallel local processes.
#
# Each batch is its own Rscript process (permulations_batch.R). Batches whose output CSV already
# exists are skipped, so an interrupted run can simply be started again.
#
# Usage (from the repo root):
#   caffeinate -i bash scripts/04_permulation_loss_dup/run_permulation_batches.sh
#
# Settings (environment variables):
#   BATCHES   number of batches (default 9; with the existing 10,000 that gives 100,000)
#   NTREES    permulations per batch (default 10000)
#   OUT       output directory (default data/perms_batches)
#   RSCRIPT   Rscript binary that has RERconverge installed (default Rscript)
#   FIRST     first batch id (default 1)
#
# Smoke test before the full run: BATCHES=1 NTREES=200 OUT=/tmp/perms_test bash <this script>

set -euo pipefail

BATCHES="${BATCHES:-9}"
NTREES="${NTREES:-10000}"
OUT="${OUT:-data/perms_batches}"
RSCRIPT="${RSCRIPT:-Rscript}"
FIRST="${FIRST:-1}"

mkdir -p "$OUT"

for ((i = FIRST; i < FIRST + BATCHES; i++)); do
    id=$(printf '%02d' "$i")
    csv="$OUT/batch_${id}_tip_values.csv"
    if [[ -s "$csv" ]]; then
        echo "Batch $id already done, skipping"
        continue
    fi
    echo "Starting batch $id ($NTREES permulations)"
    "$RSCRIPT" scripts/04_permulation_loss_dup/permulations_batch.R "$i" "$NTREES" "$OUT" \
        > "$OUT/batch_${id}.log" 2>&1 &
done

wait
echo "All batches finished. Combine with scripts/04_permulation_loss_dup/combine_permulation_batches.py"
