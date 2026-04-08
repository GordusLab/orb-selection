#!/bin/bash

# Runs BUSTED-PH with switched foreground assignment.

#SBATCH --job-name=250708_busted_ph_switch_fg
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=25,39,40,82,679,904,1077,1454,1455,1456,1809,2126,2127,2128,2129,2130
#SBATCH -n 12
#SBATCH --output=reports/%x/%A_%a.out

#make directory to store slurm reports
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REPORT_DIR=${REPORT_DIR:-$REPO_ROOT/results/slurm_reports}
mkdir -p "$REPORT_DIR/${SBATCH_JOB_NAME}/"

module load anaconda
CONDA_ENV_PATH=${CONDA_ENV_PATH:-selection}
conda activate "$CONDA_ENV_PATH"

WD=${WORK_DIR:-PATH_TO_EXTERNAL_WORK_DIR}
HOG_LIST=${HOG_LIST:-${WD}/HOG_CDS/N5.udiv.o75_list.txt}
FG_NAME=nonorb_fg	
HYPHY_ANALYSES_DIR=${HYPHY_ANALYSES_DIR:-PATH_TO_HYPHY_ANALYSES_DIR}

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

# hyphy "$HYPHY_ANALYSES_DIR/remove-duplicates/remove-duplicates.bf" \
#  CPU=3 \
#  --msa ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.fasta \
#  --tree ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.tree \
#  --output ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.no_dups.nex \
#  ENV="TOLERATE_NUMERICAL_ERRORS=1;"

# #check if BUSTED-PH has already completed for this HOG 
BUSTEDPH_OUT=${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_${FG_NAME}.json

if grep -q "p-value" ${BUSTEDPH_OUT}; then
	echo "BUSTED-PH already complete."
else
	#run busted-ph
	hyphy "$HYPHY_ANALYSES_DIR/BUSTED-PH/BUSTED-PH.bf" \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.no_dups.nex \
	 --branches FOREGROUND \
	 --output ${BUSTEDPH_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate