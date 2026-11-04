#!/bin/bash

# Runs one chunk of the phyloglm permulation null per SLURM array task
# (see phyloglm_permuted_chunk.R). Paths are placeholders -- edit before use.
#
# Usage:
#   sbatch --array=0-<n_chunks> run_phyloglm_permuted.sh OUTDIR [CHUNK_SIZE]
#     task 0   = the real labeling
#     task K>=1 = labelings (K-1)*CHUNK_SIZE+1 .. K*CHUNK_SIZE of the perms file
#   e.g. 10,000 labelings at CHUNK_SIZE=100 -> --array=0-100
#
# Resumable: finished chunks (OUTDIR/chunk_<K>.rds) are skipped, so a requeued
# array can be resubmitted as-is. Size --time/--cpus-per-task from a timing run:
# a chunk is CHUNK_SIZE x ~12,000 genes x (seconds per fit) core-seconds.
# Then run summarize_phyloglm_permuted.R --dir OUTDIR.

#SBATCH --job-name=phyloglm_permuted
#SBATCH --chdir=PATH_TO_SCRATCH
#SBATCH --partition=PATH_TO_PARTITION
#SBATCH --account=PATH_TO_ACCOUNT
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --output=PATH_TO_REPORTS/%x/%A_%a.out
#SBATCH --error=PATH_TO_REPORTS/%x/%A_%a.err

# Submit from the scratch directory with real paths on the command line (see "Launching jobs and
# storing SLURM logs" in CLAUDE.md), e.g. --chdir=$S --job-name=$(date +%y%m%d)_phyloglm_permuted
# --output=$S/reports/<job-name>/%A_%a.out --error=...%A_%a.err, after mkdir -p $S/reports/<job-name>.

set -e

OUTDIR="$1"
CHUNK_SIZE="${2:-100}"
REPO="${REPO:-PATH_TO_orb-selection}"

# R and its library path; --vanilla below skips the project .Rprofile (renv) and uses this library
module load ${R_MODULES:-PATH_TO_R_MODULES}   # e.g. R_MODULES="gfbf/2023b R/4.5.1-gfbf-2023b"
export R_LIBS_USER="${R_LIBS_USER:-PATH_TO_R_LIBRARY}"

# one BLAS thread per R worker: the chunk script already runs one worker per core, and letting
# each worker also spawn node-wide BLAS threads oversubscribes the CPUs (~8x slower measured)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$REPO"
Rscript --vanilla scripts/05_phyloglm/phyloglm_permuted_chunk.R \
	--chunk "$SLURM_ARRAY_TASK_ID" \
	--chunk-size "$CHUNK_SIZE" \
	--cores "${SLURM_CPUS_PER_TASK}" \
	--outdir "$OUTDIR" \
	${EXTRA_ARGS:-}   # e.g. EXTRA_ARGS="--perms data/perms_tip_values_500k.csv --genes FILE" to extend a subset
