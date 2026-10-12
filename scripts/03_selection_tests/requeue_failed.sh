#!/bin/bash

# Resubmits failed array tasks of a job with more CPUs and/or more time.
# Usage: requeue_failed.sh <array_job_id> <hog_list> <cpus> [time] [state_pattern]
# The original job's HOG list is needed because task IDs index into it.
# state_pattern is an awk regex matched against sacct's State column; default
# matches OUT_OF_MEMORY, TIMEOUT and FAILED. Memory is 4 GB per CPU on parallel.

set -e

JOB="$1"; LIST="$2"; CPUS="$3"; TIME="${4:-3-00:00:00}"; STATE_PATTERN="${5:-OUT_OF_ME|TIMEOUT|FAILED}"
SCRIPT="/home/crunnel2/orb-selection/scripts/03_selection_tests/run_selection_tests.sh"

IDS=$(sacct -j "$JOB" -X -n -P --format=JobID,State \
	| awk -F'|' -v pat="$STATE_PATTERN" '$2 ~ pat {split($1, a, "_"); print a[2]}' | paste -sd, -)

if [ -z "$IDS" ]; then
	echo "No tasks matching /$STATE_PATTERN/ in $JOB"
	exit 0
fi

JOB_NAME="$(date +%y%m%d)_requeue_${JOB}"
mkdir -p "/scratch4/agordus1/crunnel2/reports/${JOB_NAME}"
echo "Resubmitting tasks $IDS at $CPUS CPUs, $TIME"

sbatch --job-name="$JOB_NAME" --cpus-per-task="$CPUS" --time="$TIME" \
	--array="$IDS" --export=ALL,HOG_LIST="$LIST" "$SCRIPT"
