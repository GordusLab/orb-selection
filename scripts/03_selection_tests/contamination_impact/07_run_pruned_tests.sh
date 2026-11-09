#!/bin/bash

# One RELAX or BUSTED-PH test on a HOG with the suspected-contaminant tips pruned (one SLURM array task per HOG).
# Original HOG directories are only read, never changed. Resumable: each step is skipped if its output exists.
# Required environment (pass with sbatch --export=ALL,VAR=...; see 08_submit_pruned.sh):
#   TEST          relax | bustedph_orb | bustedph_nonorb
#   HOG_LIST      text file, one HOG per line (array task N = line N)
#   MANIFEST      prune_manifest.tsv (HOG <tab> semicolon-separated tips to remove)
#   WD_ORIG       directory with the original per-HOG folders (alignment + IQ-TREE treefile)
#   WD_OUT        directory for the pruned inputs and outputs (one folder per HOG)
#   REPO          clone containing data/orbweavers-list.txt, data/non-orbweavers-list.txt and this folder
#   CONDA_ENV     conda env with hyphy and biopython
#   HYPHY_ANALYSES  clone of veg/hyphy-analyses (LabelTrees/label-tree.bf)
# Paths below in #SBATCH lines are placeholders; override on the command line.

#SBATCH --job-name=YYMMDD_pruned_tests
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --ntasks=1
#SBATCH --chdir=PATH_TO_SCRATCH
#SBATCH --output=PATH_TO_REPORTS/%x/%A_%a.out
#SBATCH --error=PATH_TO_REPORTS/%x/%A_%a.err

set -e
: "${TEST:?}" "${HOG_LIST:?}" "${MANIFEST:?}" "${WD_ORIG:?}" "${WD_OUT:?}" "${REPO:?}" "${CONDA_ENV:?}" "${HYPHY_ANALYSES:?}"

module load anaconda
conda activate "$CONDA_ENV"
export PATH="${CONDA_ENV}/bin:${PATH}"   # 'conda activate' alone left a python without biopython on PATH in the first submission
"${CONDA_ENV}/bin/python" -c "import Bio" || { echo "biopython missing in $CONDA_ENV" >&2; exit 1; }
command -v hyphy >/dev/null || { echo "hyphy not found" >&2; exit 1; }

HOG="$(sed "${SLURM_ARRAY_TASK_ID}q;d" "$HOG_LIST")"
[ -n "$HOG" ] || { echo "no HOG for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
O="${WD_OUT}/${HOG}"; mkdir -p "$O"
ALN="${O}/${HOG}_NT.trim.dedup.pruned.nex"
TREE="${O}/${HOG}.pruned.treefile"
ORB_TREE="${O}/${HOG}.orb_fg.pruned.tree"
NONORB_TREE="${O}/${HOG}.nonorb_fg.pruned.tree"

# 1. prune alignment and tree
if [ ! -s "$ALN" ] || [ ! -s "$TREE" ]; then
	REMOVE="$(awk -F'\t' -v h="$HOG" '$1==h {print $2}' "$MANIFEST")"
	[ -n "$REMOVE" ] || { echo "no tips to remove for $HOG in $MANIFEST" >&2; exit 1; }
	"${CONDA_ENV}/bin/python" "${REPO}/scripts/03_selection_tests/contamination_impact/06_prune_hog.py" \
		--nex "${WD_ORIG}/${HOG}/${HOG}_NT.trim.dedup.nex" --tree "${WD_ORIG}/${HOG}/iqtree/${HOG}.treefile" \
		--remove "$REMOVE" --out-nex "$ALN" --out-tree "$TREE"
fi

# 2. label foreground branches exactly as in prep_for_hyphy2.sh
label () { # list tree
	local REGEX="^($(awk 'NF { sub(/\r$/, ""); printf "%s%s", sep, $0; sep = "|" }' "$1"))\\|"
	hyphy "${HYPHY_ANALYSES}/LabelTrees/label-tree.bf" --tree "$TREE" --regexp "$REGEX" --output "$2" --internal-nodes "All descendants"
}
[ -s "$ORB_TREE" ] || label "${REPO}/data/orbweavers-list.txt" "$ORB_TREE"
[ -s "$NONORB_TREE" ] || label "${REPO}/data/non-orbweavers-list.txt" "$NONORB_TREE"
for f in "$ALN" "$ORB_TREE" "$NONORB_TREE"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done

# 3. test (same options as run_selection_tests.sh)
case "$TEST" in
	bustedph_orb|bustedph_nonorb)
		if [ "$TEST" = bustedph_orb ]; then T="$ORB_TREE"; N=orb_fg; else T="$NONORB_TREE"; N=non_orb_fg; fi
		OUT="${O}/${HOG}_BUSTED-PH_${N}.pruned.json"
		if grep -qs "Comparative selection test results" "$OUT"; then echo "already complete"; else
			hyphy busted-ph "CPU=${SLURM_CPUS_PER_TASK}" --alignment "$ALN" --tree "$T" --branches Foreground \
				--multiple-hits Double+Triple --srv Yes --error-sink Yes \
				--intermediate-fits "${O}/${HOG}_BUSTED-PH_${N}.pruned_intermediate.json" --output "$OUT" ENV="TOLERATE_NUMERICAL_ERRORS=1;"
		fi ;;
	relax)
		OUT="${O}/${HOG}_RELAX.pruned.json"
		if grep -qs "relaxation or intensification parameter" "$OUT"; then echo "already complete"; else
			hyphy relax "CPU=${SLURM_CPUS_PER_TASK}" --alignment "$ALN" --tree "$NONORB_TREE" --test Foreground \
				--multiple-hits Double+Triple --srv Yes --error-sink Yes \
				--intermediate-fits "${O}/${HOG}_RELAX.pruned_intermediate.json" --output "$OUT" ENV="TOLERATE_NUMERICAL_ERRORS=1;"
		fi ;;
	*) echo "bad TEST=$TEST" >&2; exit 1 ;;
esac
conda deactivate
