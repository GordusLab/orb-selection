#!/bin/bash

# Submits 07_run_pruned_tests.sh per runtime tier (A/B/C, as in submit_tiers.sh) for RELAX and BUSTED-PH orb fg.
# Usage: 08_submit_pruned.sh [--dry-run] [--only HOG_LIST_FILE]   (tier lists: results/hyphy_contamination_impact/prune_tier_{A,B,C}.txt)
# CPUs/time follow submit_tiers.sh (RELAX 4/6/6 CPUs; BUSTED-PH 2 CPUs; 4 GB per CPU on parallel).
# Edit the variables below for the cluster; Calvin's real values were passed in the 261109 submission (see REPORT).
set -e
SCRATCH="${SCRATCH:?}"; REPO="${REPO:?}"; CONDA_ENV="${CONDA_ENV:?}"; HYPHY_ANALYSES="${HYPHY_ANALYSES:?}"
WD_ORIG="${WD_ORIG:?}"; WD_OUT="${WD_OUT:?}"
DRY=0; [ "$1" = "--dry-run" ] && DRY=1
PREFIX="$(date +%y%m%d)_pruned"
L="${REPO}/results/hyphy_contamination_impact"
for entry in "A 12:00:00 4" "B 1-12:00:00 6" "C 3-00:00:00 6"; do
	read -r TIER TIME RCPUS <<< "$entry"
	LIST="${L}/prune_tier_${TIER}.txt"; N=$(wc -l < "$LIST")
	for t in "relax $RCPUS" "bustedph_orb 2"; do
		read -r TEST CPUS <<< "$t"
		JOB="${PREFIX}_tier${TIER}_${TEST}"
		mkdir -p "${SCRATCH}/reports/${JOB}"
		CMD=(sbatch --job-name="$JOB" --chdir="$SCRATCH" --cpus-per-task="$CPUS" --time="$TIME" --array="1-${N}"
			--output="${SCRATCH}/reports/%x/%A_%a.out" --error="${SCRATCH}/reports/%x/%A_%a.err"
			--export=ALL,TEST="$TEST",HOG_LIST="$LIST",MANIFEST="${L}/prune_manifest.tsv",WD_ORIG="$WD_ORIG",WD_OUT="$WD_OUT",REPO="$REPO",CONDA_ENV="$CONDA_ENV",HYPHY_ANALYSES="$HYPHY_ANALYSES"
			"${REPO}/scripts/03_selection_tests/contamination_impact/07_run_pruned_tests.sh")
		if [ $DRY = 1 ]; then echo "${CMD[@]}"; else "${CMD[@]}"; fi
	done
done
