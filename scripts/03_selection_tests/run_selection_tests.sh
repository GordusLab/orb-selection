#!/bin/bash

# Runs RELAX and BUSTED-PH (orb fg and non orb fg) for one HOG per SLURM array task.

#SBATCH --array=1-4756
#SBATCH -n 3
#SBATCH --output=reports/%x/%A_%a.out

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/hyphy-new

WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt

#HYPHY_ANALYSES_DIR=/home/crunnel2/bin/hyphy-analyses/

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

## RELAX
RELAX_OUT=${WD}/${CURRENT_HOG}_RELAX2.json
#check if RELAX has already completed for this HOG 
if grep -q "p-value" ${RELAX_OUT}; then
	echo "RELAX already complete."
else
	#run relax
	hyphy relax \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.non_orb_fg.nex \
	 --test "FOREGROUND" \
	 --multiple-hits Double+Triple \
	 --models Minimal \
	 --srv Yes \
 	 ENV="TOLERATE_NUMERICAL_ERRORS=1;" \
	 --output ${RELAX_OUT}
fi

## BUSTED-PH

BUSTEDPH_ORB_OUT=${WD}/${CURRENT_HOG}_BUSTED-PH_orb_fg.json
# check if BUSTED-PH, orb fg has already completed for this HOG 
if grep -q "p-value" ${BUSTEDPH_ORB_OUT}; then
	echo "BUSTED-PH-fw already complete."
else
	#run busted-ph
	hyphy busted-ph \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.orb_fg.nex \
	 --branches FOREGROUND \
	 --output ${BUSTEDPH_ORB_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

BUSTEDPH_NON_ORB_OUT=${WD}/${CURRENT_HOG}_BUSTED-PH_non_orb_fg.json
# check if BUSTED-PH, non-orb fg has already completed for this HOG 
if grep -q "p-value" ${BUSTEDPH_NON_ORB_OUT}; then
	echo "BUSTED-PH-rev already complete."
else
	#run busted-ph
	hyphy busted-ph \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.non_orb_fg.nex \
	 --branches FOREGROUND \
	 --output ${BUSTEDPH_NON_ORB_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate
