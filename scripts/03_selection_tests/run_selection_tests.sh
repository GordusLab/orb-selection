#!/bin/bash

# Runs one RELAX or BUSTED-PH test per SLURM array task, with three tasks per HOG.

#SBATCH --job-name=261002_run_selection_tests
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=06:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-14268
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --output=/scratch4/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/scratch4/agordus1/crunnel2/reports/%x/%A_%a.err

set -e

# Make directory to store SLURM reports
mkdir -p "/scratch4/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/"

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/hyphy-new

WD="/scratch4/agordus1/crunnel2/hyphy_wd_260929"
HOG_LIST="${HOG_LIST:-/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt}"

HOG_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / 3 + 1 ))
TEST_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % 3 + 1 ))
CURRENT_HOG="$(sed "${HOG_INDEX}q;d" "$HOG_LIST")"

ALN_FILE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_NT.trim.dedup.nex"
ORB_TREE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}.orb_fg.tree"
NONORB_TREE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}.nonorb_fg.tree"

for f in "$ALN_FILE" "$ORB_TREE" "$NONORB_TREE"; do
	if [ ! -s "$f" ]; then
		echo "ERROR: missing or empty input for ${CURRENT_HOG}: $f" >&2
		exit 1
	fi
done

case "$TEST_INDEX" in
	1)
		BUSTEDPH_ORB_OUT="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_orb_fg.json"
		if grep -q "p-value" "$BUSTEDPH_ORB_OUT"; then
			echo "BUSTED-PH, orb fg already complete."
		else
			hyphy busted-ph \
				"CPU=${SLURM_CPUS_PER_TASK}" \
				--alignment "$ALN_FILE" \
				--tree "$ORB_TREE" \
				--branches Foreground \
				--multiple-hits Double+Triple \
				--srv Yes \
				--error-sink Yes \
				--intermediate-fits "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_orb_fg_intermediate.json" \
				--output "$BUSTEDPH_ORB_OUT" \
				ENV="TOLERATE_NUMERICAL_ERRORS=1;"
		fi
		;;
	2)
		BUSTEDPH_NON_ORB_OUT="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_non_orb_fg.json"
		if grep -q "p-value" "$BUSTEDPH_NON_ORB_OUT"; then
			echo "BUSTED-PH, non-orb fg already complete."
		else
			hyphy busted-ph \
				"CPU=${SLURM_CPUS_PER_TASK}" \
				--alignment "$ALN_FILE" \
				--tree "$NONORB_TREE" \
				--branches Foreground \
				--multiple-hits Double+Triple \
				--srv Yes \
				--error-sink Yes \
				--intermediate-fits "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH_non_orb_fg_intermediate.json" \
				--output "$BUSTEDPH_NON_ORB_OUT" \
				ENV="TOLERATE_NUMERICAL_ERRORS=1;"
		fi
		;;
	3)
		RELAX_OUT="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX.json"
		if grep -q "p-value" "$RELAX_OUT"; then
			echo "RELAX already complete."
		else
			hyphy relax \
				"CPU=${SLURM_CPUS_PER_TASK}" \
				--alignment "$ALN_FILE" \
				--tree "$NONORB_TREE" \
				--test Foreground \
				--multiple-hits Double+Triple \
				--srv Yes \
				--error-sink Yes \
				--intermediate-fits "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX_intermediate.json" \
				--output "$RELAX_OUT" \
				ENV="TOLERATE_NUMERICAL_ERRORS=1;"
		fi
		;;
	esac

conda deactivate
