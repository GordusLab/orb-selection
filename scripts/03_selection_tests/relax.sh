#!/bin/bash

# Runs RELAX for one HOG per SLURM array task.

#SBATCH --array=1-4756
#SBATCH -n 12
#SBATCH --output=reports/%x/%A_%a.out

#make directory to store slurm reports
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REPORT_DIR=${REPORT_DIR:-$REPO_ROOT/results/slurm_reports}
mkdir -p "$REPORT_DIR/$SBATCH_JOB_NAME/"

module load anaconda
CONDA_ENV_PATH=${CONDA_ENV_PATH:-selection}
conda activate "$CONDA_ENV_PATH"

WD=${WORK_DIR:-PATH_TO_EXTERNAL_WORK_DIR}
HOG_LIST=${HOG_LIST:-${WD}/HOG_CDS/N5.udiv.o75_list.txt}
FG_NAME=ow_fg_cons

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
	CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "$HOG_LIST")
else
	CURRENT_HOG=${CURRENT_HOG:-PATH_TO_HOG_ID}
fi

#check if RELAX has already completed for this HOG 
RELAX_OUT=${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX.json

if grep -q "p-value" ${RELAX_OUT}; then
	echo "RELAX already complete."
else
	#run relax
	hyphy relax \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.fasta \
	 --tree ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.tree \
	 --test "Unlabeled branches" \
	 --multiple-hits Double+Triple \
	 --models Minimal \
	 --srv Yes \
 	 ENV="TOLERATE_NUMERICAL_ERRORS=1;" \
	 --output ${RELAX_OUT}
fi

conda deactivate
