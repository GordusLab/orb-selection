#!/bin/bash

# Runs BUSTED-PH with switched foreground assignment.

#SBATCH --array=1-4756
#SBATCH -n 12
#SBATCH --output=reports/%x/%A_%a.out

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

#make directory to store slurm reports
REPORT_DIR=PATH_TO_REPORT_DIR
mkdir -p "$REPORT_DIR/${SBATCH_JOB_NAME}/"

module load anaconda
CONDA_ENV_PATH=PATH_TO_CONDA_ENV
conda activate "$CONDA_ENV_PATH"

WD=PATH_TO_EXTERNAL_WORK_DIR
HOG_LIST=${HOG_LIST:-$REPO_ROOT/data/N5.udiv.o75_list.txt}
FG_NAME=nonorb_fg	

#path to hyphy analyses repo clone: https://github.com/veg/hyphy-analyses/
HYPHY_ANALYSES_DIR=PATH_TO_HYPHY_ANALYSES_DIR_CLONE

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

# #check if BUSTED-PH has already completed for this HOG 
BUSTEDPH_OUT=${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_${FG_NAME}.json

if grep -q "p-value" ${BUSTEDPH_OUT}; then
	echo "BUSTED-PH already complete."
else
	#run busted-ph
	hyphy "$HYPHY_ANALYSES_DIR/BUSTED-PH/BUSTED-PH.bf" \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.no_dups.nex \
	 --branches "Unlabeled branches" \
	 --output ${BUSTEDPH_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate