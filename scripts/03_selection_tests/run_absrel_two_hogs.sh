#!/bin/bash

# Runs aBSREL for two specific HOGs as a SLURM array job.

#SBATCH --job-name=261020_absrel_non_orb
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=24:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=NONE
#SBATCH --array=1-21
#SBATCH --chdir=/scratch4/agordus1/crunnel2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --output=/scratch4/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/scratch4/agordus1/crunnel2/reports/%x/%A_%a.err

mkdir -p /scratch4/agordus1/crunnel2/reports/261020_absrel_non_orb

set -e

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/hyphy-new

WD="/scratch4/agordus1/crunnel2/hyphy_wd_260929"
HOG_LIST="/scratch4/agordus1/crunnel2/nonorb_fg.txt"

HOG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$HOG_LIST")

echo "Running aBSREL for $HOG (task $SLURM_ARRAY_TASK_ID)..."

ALN_FILE="${WD}/${HOG}/${HOG}_NT.trim.dedup.nex"
NONORB_TREE="${WD}/${HOG}/${HOG}.nonorb_fg.tree"
ABSREL_OUT="${WD}/${HOG}/${HOG}_ABSREL_non_orb_fg.json"

# Check inputs exist
for f in "$ALN_FILE" "$NONORB_TREE"; do
	if [ ! -s "$f" ]; then
		echo "ERROR: missing or empty input for ${HOG}: $f" >&2
		exit 1
	fi
done

# Skip if already complete
if grep -qs "p-value\|overall" "$ABSREL_OUT" 2>/dev/null; then
	echo "aBSREL for $HOG already complete."
else
	hyphy absrel \
		"CPU=${SLURM_CPUS_PER_TASK}" \
		--alignment "$ALN_FILE" \
		--tree "$NONORB_TREE" \
		--multiple-hits Double+Triple \
		--srv Yes \
		--error-sink Yes \
		--intermediate-fits "${WD}/${HOG}/${HOG}_ABSREL_non_orb_fg_intermediate.json" \
		--output "$ABSREL_OUT" \
		ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

echo "Done: $ABSREL_OUT"

conda deactivate
