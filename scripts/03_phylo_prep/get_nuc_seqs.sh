#!/bin/bash

# Extract coding nucleotide sequences for one HOG per array task.

#SBATCH --job-name=260304_get_nuc_seqs
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=00:05:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-3375
#SBATCH --output=reports/%x/%A_%a.out

#run the program from desired CDS output directory

module load anaconda
CONDA_ENV_PATH=${CONDA_ENV_PATH:-seqkit}
conda activate "$CONDA_ENV_PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

#make directory to store slurm reports
REPORT_DIR=${REPORT_DIR:-$REPO_ROOT/results/slurm_reports}
mkdir -p "$REPORT_DIR/$SBATCH_JOB_NAME/"

HOG_LIST=${HOG_LIST:-$REPO_ROOT/assets/N5_occ48_remaining.txt}
NX_TSV=${NX_TSV:-$REPO_ROOT/data/N5.tsv}
TD_DIR=${TD_DIR:-PATH_TO_TRANSDECODER_DIR}

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

python "$SCRIPT_DIR/get_nuc_seqs.py" "$NX_TSV" "$CURRENT_HOG" "$TD_DIR"

conda deactivate