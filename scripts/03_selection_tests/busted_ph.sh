#!/bin/bash

# Runs BUSTED-PH for one HOG per SLURM array task.
# Argument: foreground assignment name (e.g., "orb_fg" or "nonorb_fg")

#SBATCH --job-name=260722_busted_ph_rev_test
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=00:10:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=4746-4756
#SBATCH -n 12
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt
FG_NAME=$1

HYPHY_ANALYSES_DIR=/home/crunnel2/bin/hyphy-analyses/

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

# check if BUSTED-PH has already completed for this HOG 
BUSTEDPH_OUT=${WD}/${CURRENT_HOG}_BUSTED-PH_${FG_NAME}.json

if [ $FG_NAME == "orb_fg" ]; then
	BRANCHES="Foreground"
else
	BRANCHES="Unlabeled branches"
fi

if grep -q "p-value" ${BUSTEDPH_OUT}; then
	echo "BUSTED-PH already complete."
else
	#run busted-ph
	hyphy "$HYPHY_ANALYSES_DIR/BUSTED-PH/BUSTED-PH.bf" \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.fltrd.fasta \
	 --tree ${WD}/${CURRENT_HOG}.orb_fg.tree \
	 --branches ${BRANCHES} \
	 --output ${BUSTEDPH_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate