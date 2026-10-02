#!/bin/bash

# Runs RELAX and BUSTED-PH (orb fg and non-orb fg) for one HOG per SLURM array task.

#SBATCH --job-name=261002_run_selection_tests_hyphymp
#SBATCH --partition=shared
#SBATCH --account=agordus1
#SBATCH --time=06:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=3-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

set -e

# Make directory to store SLURM reports
mkdir -p "/data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/"

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/hyphy-new

WD="/scratch4/agordus1/crunnel2/hyphy_wd_260929"
HOG_LIST="/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt"

#HYPHY_ANALYSES_DIR=/home/crunnel2/bin/hyphy-analyses/

CURRENT_HOG="$(sed "${SLURM_ARRAY_TASK_ID}q;d" "$HOG_LIST")"

ALN_FILE="${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.dedup.nex"
ORB_TREE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}.orb_fg.tree"
NONORB_TREE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}.nonorb_fg.tree"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

## BUSTED-PH

BUSTEDPH_ORB_OUT="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_orb_fg.json"
# Check if BUSTED-PH, orb fg has already completed for this HOG
if grep -q "p-value" "$BUSTEDPH_ORB_OUT"; then
	echo "BUSTED-PH, orb fg already complete."
else
	# Run BUSTED-PH
	hyphy busted-ph \
		--alignment "$ALN_FILE" \
		--tree "$ORB_TREE" \
		--branches Foreground \
		--multiple-hits Double+Triple \
		--srv Yes \
		--error-sink Yes \
		--intermediate-fits "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_orb_fg_intermediate.json" \
		--output "$BUSTEDPH_ORB_OUT" \
		CPU=${SLURM_CPUS_PER_TASK} \
		ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

BUSTEDPH_NON_ORB_OUT="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_non_orb_fg.json"
# Check if BUSTED-PH, non-orb fg has already completed for this HOG
if grep -q "p-value" "$BUSTEDPH_NON_ORB_OUT"; then
	echo "BUSTED-PH, non-orb fg already complete."
else
	# Run BUSTED-PH
	hyphy busted-ph \
		--alignment "$ALN_FILE" \
		--tree "$NONORB_TREE" \
		--branches Foreground \
		--multiple-hits Double+Triple \
		--srv Yes \
		--error-sink Yes \
		--intermediate-fits "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_non_orb_fg_intermediate.json" \
		--output "$BUSTEDPH_NON_ORB_OUT" \
		CPU=${SLURM_CPUS_PER_TASK} \
		ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

## RELAX
RELAX_OUT="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX.json"
# Check if RELAX has already completed for this HOG
if grep -q "p-value" "$RELAX_OUT"; then
	echo "RELAX already complete."
else
	# Run RELAX
	hyphy relax \
		--alignment "$ALN_FILE" \
		--tree "$NONORB_TREE" \
		--test Foreground \
		--multiple-hits Double+Triple \
		--srv Yes \
		--error-sink Yes \
		--intermediate-fits "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX_intermediate.json" \
		--output "$RELAX_OUT" \
		CPU=${SLURM_CPUS_PER_TASK} \
		ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate
