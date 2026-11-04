#!/bin/bash

# Runs RELAX on one (HOG, permuted labeling) pair per SLURM array task, as a null
# calibration for the real RELAX results (see make_permuted_relax_manifest.py).
# Reuses each HOG's existing alignment and IQ-TREE tree from the real run; only the
# foreground labeling changes. Arguments/paths are placeholders -- edit before use.
#
# Usage:
#   sbatch --array=1-<n_tasks> run_relax_permuted.sh MANIFEST WD OUT_ROOT HYPHY_ANALYSES
#     MANIFEST        manifest.tsv from make_permuted_relax_manifest.py (HOG, perm_id, regex)
#     WD              work dir holding <HOG>/<HOG>_NT.trim.dedup.nex and <HOG>/iqtree/<HOG>.treefile
#     OUT_ROOT        output root; writes <OUT_ROOT>/<perm_id>/<HOG>_RELAX.json
#     HYPHY_ANALYSES  hyphy-analyses checkout (for LabelTrees/label-tree.bf)
#
# Set CONDA_ENV to the HyPhy conda env (e.g. via sbatch --export=ALL,CONDA_ENV=...), and override
# --partition/--account/--output/--time/--cpus-per-task on the sbatch command line as needed.
#
# Resumable like the stage 02/03 scripts: finished outputs (JSON containing
# "p-value") are skipped, so a requeued task picks up where it left off. Size
# --time/--cpus-per-task like the real RELAX tier for the sampled HOGs
# (see submit_tiers.sh).

#SBATCH --job-name=relax_permuted
#SBATCH --chdir=PATH_TO_SCRATCH
#SBATCH --partition=PATH_TO_PARTITION
#SBATCH --account=PATH_TO_ACCOUNT
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=PATH_TO_REPORTS/%x/%A_%a.out
#SBATCH --error=PATH_TO_REPORTS/%x/%A_%a.err

# Submit from the scratch directory with real paths on the command line (see "Launching jobs and
# storing SLURM logs" in CLAUDE.md), e.g. --chdir=$S --job-name=$(date +%y%m%d)_relax_permuted
# --output=$S/reports/<job-name>/%A_%a.out --error=...%A_%a.err, after mkdir -p $S/reports/<job-name>.

set -e

MANIFEST="$1"
WD="$2"
OUT_ROOT="$3"
HYPHY_ANALYSES="$4"

module load anaconda
conda activate "${CONDA_ENV:-PATH_TO_HYPHY_CONDA_ENV}"

IFS=$'\t' read -r CURRENT_HOG PERM_ID REGEX < <(sed "${SLURM_ARRAY_TASK_ID}q;d" "$MANIFEST")

ALN_FILE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_NT.trim.dedup.nex"
IQTREE_FILE="${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}.treefile"
for f in "$ALN_FILE" "$IQTREE_FILE"; do
	if [ ! -s "$f" ]; then
		echo "ERROR: missing or empty input for ${CURRENT_HOG}: $f" >&2
		exit 1
	fi
done

PERM_DIR="${OUT_ROOT}/${PERM_ID}"
mkdir -p "$PERM_DIR"
PERM_TREE="${PERM_DIR}/${CURRENT_HOG}.nonorb_fg.tree"
RELAX_OUT="${PERM_DIR}/${CURRENT_HOG}_RELAX.json"

if grep -qs "p-value" "$RELAX_OUT"; then
	echo "RELAX for ${CURRENT_HOG} under ${PERM_ID} already complete."
	exit 0
fi

# label the permuted foreground (same label-tree.bf call as prep_for_hyphy2.sh)
if [ ! -s "$PERM_TREE" ]; then
	hyphy "${HYPHY_ANALYSES}/LabelTrees/label-tree.bf" \
		--tree "$IQTREE_FILE" \
		--regexp "$REGEX" \
		--output "$PERM_TREE" \
		--internal-nodes "All descendants"
fi

# a regex that matched nothing would silently run RELAX with no foreground
N_FG=$(grep -o "{Foreground}" "$PERM_TREE" | wc -l)
if [ "$N_FG" -eq 0 ]; then
	echo "ERROR: no Foreground branches labeled in ${PERM_TREE}" >&2
	exit 1
fi
echo "${CURRENT_HOG} ${PERM_ID}: ${N_FG} Foreground-labeled nodes"

hyphy relax \
	"CPU=${SLURM_CPUS_PER_TASK}" \
	--alignment "$ALN_FILE" \
	--tree "$PERM_TREE" \
	--test Foreground \
	--multiple-hits Double+Triple \
	--srv Yes \
	--error-sink Yes \
	--intermediate-fits "${PERM_DIR}/${CURRENT_HOG}_RELAX_intermediate.json" \
	--output "$RELAX_OUT" \
	ENV="TOLERATE_NUMERICAL_ERRORS=1;"

conda deactivate
