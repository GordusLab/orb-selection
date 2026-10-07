#!/bin/bash

# Resubmits OUT_OF_MEMORY array tasks of a job with more CPUs (4 GB each on parallel).
# Usage: requeue_oom.sh <array_job_id> <hog_list> <cpus> [time]
# The original job's HOG list is needed because task IDs index into it.

set -e

JOB="$1"; LIST="$2"; CPUS="$3"; TIME="${4:-3-00:00:00}"
SCRIPT="/home/crunnel2/orb-selection/scripts/03_selection_tests/run_selection_tests.sh"

IDS=$(sacct -j "$JOB" -X -n -P --format=JobID,State \
	| awk -F'|' '$2 ~ /OUT_OF_ME/ {split($1, a, "_"); print a[2]}' | paste -sd, -)

if [ -z "$IDS" ]; then
	echo "No OUT_OF_MEMORY tasks in $JOB"
	exit 0
fi

JOB_NAME="$(date +%y%m%d)_requeue_oom_${JOB}"
mkdir -p "/scratch4/agordus1/crunnel2/reports/${JOB_NAME}"
echo "Resubmitting tasks $IDS at $CPUS CPUs"

sbatch --job-name="$JOB_NAME" --cpus-per-task="$CPUS" --time="$TIME" \
	--array="$IDS" --export=ALL,HOG_LIST="$LIST" "$SCRIPT"
