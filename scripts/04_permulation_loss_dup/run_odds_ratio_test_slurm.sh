#!/bin/bash
#SBATCH --job-name=odds_ratio_test
#SBATCH --chdir=PATH_TO_SCRATCH
#SBATCH --partition=shared
#SBATCH --account=agordus1
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=PATH_TO_REPORTS/%x/%j.out
#SBATCH --error=PATH_TO_REPORTS/%x/%j.err

# Run the odds ratio permulation test (run_odds_ratio_test.py) as a single SLURM job.
# With the batched implementation, 100,000 permulations take roughly 10 minutes on one core;
# the 32 GB covers the permulated tip assignments held in memory and the pickled results.
#
# Usage: submit from the scratch directory, with the real paths on the command line (the placeholders
# above are overridden; see "Launching jobs and storing SLURM logs" in CLAUDE.md). Anything after the
# script name is passed on to run_odds_ratio_test.py (see its --help):
#   S=/scratch4/agordus1/crunnel2; N=$(date +%y%m%d)_odds_ratio_test; mkdir -p $S/reports/$N
#   cd $S && sbatch --chdir=$S --job-name=$N --output=$S/reports/$N/%j.out --error=$S/reports/$N/%j.err \
#       --export=ALL,REPO_DIR=$HOME/orb-selection \
#       $HOME/orb-selection/scripts/04_permulation_loss_dup/run_odds_ratio_test_slurm.sh \
#       --permulations-csv data/perms_tip_values_100k.csv.gz --dir-suffix 100k
#
# Settings (environment variables, pass with sbatch --export=ALL,NAME=value):
#   PYTHON_BIN      python with numpy, pandas, scipy, matplotlib, seaborn and tqdm (default: python)
#   PYTHON_MODULES  modules to load first, if the python comes from a module (optional)
#   REPO_DIR        repo root (default: PATH_TO_orb-selection; the job starts in the scratch directory)

set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"
REPO_DIR="${REPO_DIR:-PATH_TO_orb-selection}"

cd "$REPO_DIR"

if [[ -n "${PYTHON_MODULES:-}" ]]; then
    # shellcheck disable=SC2086  # intentional word splitting: one module name per word
    module load $PYTHON_MODULES
fi

echo "Started $(date) with $("$PYTHON_BIN" --version 2>&1)"
"$PYTHON_BIN" scripts/04_permulation_loss_dup/run_odds_ratio_test.py "$@"
echo "Finished $(date)"
