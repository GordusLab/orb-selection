#!/bin/bash

# Runs RELAX and BUSTED-PH (fw and rev) for one HOG per SLURM array task.

#SBATCH --job-name=260804_run_selection_tests_timeout_threaded_6cpus
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=21,33,35,43,48,82,162,266,357,377,482,656,708,800,828,860,909,953,958,1011,1104,1190,1316,1533,1634,2106,2460,2493,2495,2530,2532,2542,2558,2674,2723,2823,2825,2849,2945,2988,2997,3042,3049,3075,3105,3142,3191,3365,3382,3397,3424,3447,3478,3822,3837,3838,3839,3840,3855,3871,3936,3947,4130,4286,4441,4500,4519,4529
#SBATCH -n 6
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt

HYPHY_ANALYSES_DIR=/home/crunnel2/bin/hyphy-analyses/

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

## RELAX
RELAX_OUT=${WD}/${CURRENT_HOG}_RELAX.json
#check if RELAX has already completed for this HOG 
if grep -q "p-value" ${RELAX_OUT}; then
	echo "RELAX already complete."
else
	#run relax
	hyphy relax \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.orb_fg.nex \
	 --test "Unlabeled branches" \
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
	hyphy "$HYPHY_ANALYSES_DIR/BUSTED-PH/BUSTED-PH.bf" \
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
	hyphy "$HYPHY_ANALYSES_DIR/BUSTED-PH/BUSTED-PH.bf" \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.non_orb_fg.nex \
	 --branches FOREGROUND \
	 --output ${BUSTEDPH_NON_ORB_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate
