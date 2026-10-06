#!/bin/bash

# Submits run_selection_tests.sh once per runtime tier (lists in data/hyphy_tiers/),
# overriding CPUs, time limit, array size, and HOG list. Run from the repo root.
# Tiers come from Pass 2 run times (3 CPUs, all three tests in series):
#   A < 12 h, B 12-36 h, C > 36 h / timed out / no record / known-heavy.

set -e

REPORT_ROOT="/data/agordus1/crunnel2/reports"
REPO="/home/crunnel2/orb-selection"
SCRIPT="${REPO}/scripts/03_selection_tests/run_selection_tests.sh"
PREFIX="$(date +%y%m%d)_run_selection_tests"

# tier  cpus  time
TIERS=(
	"A 3 12:00:00"
	"B 4 1-12:00:00"
	"C 6 3-00:00:00"
)

for entry in "${TIERS[@]}"; do
	read -r TIER CPUS TIME <<< "$entry"
	LIST="${REPO}/data/hyphy_tiers/tier_${TIER}.txt"
	N_TASKS=$(( 3 * $(wc -l < "$LIST") ))
	JOB_NAME="${PREFIX}_tier${TIER}"

	# SLURM opens the output files before the script runs, so the directory must exist now
	mkdir -p "${REPORT_ROOT}/${JOB_NAME}"

	sbatch \
		--job-name="$JOB_NAME" \
		--cpus-per-task="$CPUS" \
		--time="$TIME" \
		--array="1-${N_TASKS}" \
		--export=ALL,HOG_LIST="$LIST" \
		"$SCRIPT"
done
