#!/bin/bash

# Submits run_selection_tests.sh per runtime tier (lists in data/hyphy_tiers/) and per test type.
# Array task N maps to test (N-1)%3+1: 1 = BUSTED-PH orb, 2 = BUSTED-PH non-orb, 3 = RELAX,
# so each test is a strided array (start:3 step).
# Tiers come from Pass 2 run times: A < 12 h, B 12-36 h, C > 36 h / no record / known-heavy.
# Memory is capped at 4 GB per CPU on parallel, so CPUs are sized for memory:
# BUSTED-PH 2 CPUs (8 GB), RELAX 4 CPUs (16 GB). Requeue OOM tasks with requeue_oom.sh.

set -e

REPORT_ROOT="/scratch4/agordus1/crunnel2/reports"
REPO="/home/crunnel2/orb-selection"
SCRIPT="${REPO}/scripts/03_selection_tests/run_selection_tests.sh"
PREFIX="$(date +%y%m%d)_run_selection_tests"

# tier  time
TIERS=(
	"A 12:00:00"
	"B 1-12:00:00"
	"C 3-00:00:00"
)

# label  first_task_offset  cpus
TESTS=(
	"bustedph_orb 1 2"
	"bustedph_nonorb 2 2"
	"relax 3 4"
)

for entry in "${TIERS[@]}"; do
	read -r TIER TIME <<< "$entry"
	LIST="${REPO}/data/hyphy_tiers/tier_${TIER}.txt"
	N_TASKS=$(( 3 * $(wc -l < "$LIST") ))

	for t in "${TESTS[@]}"; do
		read -r LABEL START CPUS <<< "$t"
		JOB_NAME="${PREFIX}_tier${TIER}_${LABEL}"

		# SLURM opens the output files before the script runs, so the directory must exist now
		mkdir -p "${REPORT_ROOT}/${JOB_NAME}"

		sbatch \
			--job-name="$JOB_NAME" \
			--cpus-per-task="$CPUS" \
			--time="$TIME" \
			--array="${START}-${N_TASKS}:3" \
			--export=ALL,HOG_LIST="$LIST" \
			"$SCRIPT"
	done
done
