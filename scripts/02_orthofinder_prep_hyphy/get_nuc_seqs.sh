#!/bin/bash

# Extract coding nucleotide sequences for one HOG per array task.

#SBATCH --array=1-4756
#SBATCH --output=reports/%x/%A_%a.out

# Usage: scripts/02_orthofinder_prep_hyphy/get_nuc_seqs.sh <conda_env_path> <report_dir> <td_dir>
# Arguments:
#   1: path to conda environment
#   2: directory to store slurm reports
#   3: directory containing TransDecoder nucleotide outputs
# Run the program from desired HOG_CDS output directory

module load anaconda
CONDA_ENV_PATH=$1
conda activate "$CONDA_ENV_PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

#make directory to store slurm reports
REPORT_DIR=$2
mkdir -p "$REPORT_DIR/$SBATCH_JOB_NAME/"

HOG_LIST=${HOG_LIST:-$REPO_ROOT/assets/N5.udiv.o75_list.txt}
NX_TSV=${NX_TSV:-$REPO_ROOT/data/N5.tsv}
TD_DIR=$3

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

python "$SCRIPT_DIR/get_nuc_seqs.py" "$NX_TSV" "$CURRENT_HOG" "$TD_DIR"

conda deactivate