#!/bin/bash
#SBATCH --job-name=permulation_batches
#SBATCH --partition=shared
#SBATCH --account=agordus1
#SBATCH --time=01:30:00
#SBATCH --array=1-45%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M  # one CPU's share on shared (4000 MB/CPU); measured peak is ~360 MB. Larger requests are charged extra CPUs.
#SBATCH --output=permulation_batches_%A_%a.out
#SBATCH --error=permulation_batches_%A_%a.err

# One categorical-permulation batch per array task (permulations_batch.R).
# Each task is independent and resumable: a task whose output CSV already exists exits immediately,
# so a failed or requeued task can simply be resubmitted.
#
# Usage (from the repo root on the cluster):
#   sbatch scripts/04_permulation_loss_dup/permulations_batch_slurm.sh
#   # then combine locally or on the cluster:
#   python scripts/04_permulation_loss_dup/combine_permulation_batches.py --include data/perms_tip_values.csv
#
# Settings (environment variables, pass with sbatch --export=ALL,NAME=value):
#   NTREES       permulations per batch (default 2000; about 1 s per permulation on Rockfish)
#   OUT          output directory (default data/perms_batches, relative to the repo root)
#   REPO_DIR     repo root (default: the directory sbatch was run from)
#   R_MODULES    modules to load for an R that has RERconverge (default: Rockfish's gfbf/2023b
#                toolchain plus R/4.5.1-gfbf-2023b; the toolchain must be loaded before R)
#   R_LIBS_USER  user R library containing RERconverge (default: $HOME/R/x86_64-pc-linux-gnu-library/4.5)
#
# The array size is the number of batches: 45 batches of 2,000 plus the existing 10,000 = 100,000.
# At most 20 tasks run at once (the %20 throttle).

set -euo pipefail

NTREES="${NTREES:-2000}"
OUT="${OUT:-data/perms_batches}"
REPO_DIR="${REPO_DIR:-$SLURM_SUBMIT_DIR}"

cd "$REPO_DIR"

R_MODULES="${R_MODULES:-gfbf/2023b R/4.5.1-gfbf-2023b}"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.5}"

# shellcheck disable=SC2086  # intentional word splitting: one module name per word
module load $R_MODULES

id=$(printf '%02d' "$SLURM_ARRAY_TASK_ID")
csv="$OUT/batch_${id}_tip_values.csv"
mkdir -p "$OUT"

if [[ -s "$csv" ]]; then
    echo "Batch $id already done, skipping"
    exit 0
fi

echo "Batch $id: $NTREES permulations, started $(date)"
# --vanilla skips the repo's .Rprofile, which would activate renv and hide the user library holding RERconverge
Rscript --vanilla scripts/04_permulation_loss_dup/permulations_batch.R "$SLURM_ARRAY_TASK_ID" "$NTREES" "$OUT"
echo "Batch $id finished $(date)"
