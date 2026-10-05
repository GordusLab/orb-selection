#!/bin/bash

# Prepares HOG alignments/trees for HyPhy analyses.
# Usage: scripts/02_orthofinder_prep_hyphy/prep_for_hyphy2.sh <conda_env_path> <work_dir> <hog_cds_dir> <macse_jar> <hyphy_analyses>
# Note: The script assumes that the input HOG CDS files are in FASTA format and that MACSE is available at the specified JAR path.
# Assumes that the repo is located at $HOME/orb-selection

# sbatch ~/orb-selection/scripts/02_orthofinder_prep_hyphy/prep_for_hyphy2.sh ~/anaconda3/envs/hyphy-new/ ~/scratch/hyphy_wd_260929 ~/scratch/hyphy_wd_260929/HOG_CDS/ ~/bin/macse_v2.07.jar ~/bin/hyphy-analyses/

#SBATCH --job-name=261005_prep_for_hyphy2_clipkit_test
#SBATCH --partition=shared
#SBATCH --account=agordus1
#SBATCH --time=00:01:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1,2,21,82,630,1053
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

set -e

# Make directory to store SLURM reports
mkdir -p "/data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/"

CONDA_ENV_PATH="$1"
WD="$2"
HOG_CDS_DIR="$3"
MACSE_JAR="$4"
HYPHY_ANALYSES="$5"

module load anaconda
conda activate "$CONDA_ENV_PATH"

REPO_ROOT="$HOME/orb-selection"

# List of HOG IDs from directory of FASTA files
HOG_LIST="${REPO_ROOT}/data/N5.udiv.o75_list.txt"

CURRENT_HOG="$(sed "${SLURM_ARRAY_TASK_ID}q;d" "$HOG_LIST")"
mkdir -p "${WD}/${CURRENT_HOG}/"

# CDS file needs to be in prequal folder so prequal will output to that folder
CDS_FILE="${WD}/${CURRENT_HOG}/prequal/${CURRENT_HOG}.fasta"

# Check if the CDS file has already been moved to the prequal directory
if [ -f "$CDS_FILE" ]; then
	echo "CDS file ${CDS_FILE} already moved; on to prequal."
else
	mkdir -p "${WD}/${CURRENT_HOG}/prequal/"
	mv "${HOG_CDS_DIR}/${CURRENT_HOG}.fasta" "${WD}/${CURRENT_HOG}/prequal/"
fi

#############
## PREQUAL ##
#############

# Check if prequal has already completed for this HOG
PREQUAL_FILE="${WD}/${CURRENT_HOG}/prequal/${CURRENT_HOG}.fasta.dna.filtered"

if [ -f "$PREQUAL_FILE" ]; then
	echo "Prequal file ${PREQUAL_FILE} exists; on to macse."
else
	# Run prequal
	prequal -dosummary "$CDS_FILE"
fi

###########
## MACSE ##
###########

# Check if MACSE has already completed for this HOG
MACSE_NT_FILE="${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.fasta"

if [ -f "$MACSE_NT_FILE" ]; then
	echo "Macse file ${MACSE_NT_FILE} exists; on to ClipKIT."
else
	# Run MACSE
	mkdir -p "${WD}/${CURRENT_HOG}/macse/"
	java -jar -Xmx47G "$MACSE_JAR" \
		-prog alignSequences \
		-seq "$PREQUAL_FILE" \
		-out_NT "$MACSE_NT_FILE" \
		-out_AA "${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_AA.fasta"
fi

##############
## CLIPKIT  ##
##############

# Trim the MACSE codon alignment with ClipKIT. Start with smart-gap; if it picks a
# gaps threshold > 0.9 (too permissive), rerun with gappy -g 0.9.
CLIPKIT_DIR="${WD}/${CURRENT_HOG}/clipkit"
TRIM_FILE="${CLIPKIT_DIR}/${CURRENT_HOG}_NT.trim.fasta"
CLIPKIT_STATS="${CLIPKIT_DIR}/${CURRENT_HOG}_clipkit_stats.tsv"

# Set RERUN_FROM_TRIM=1 to redo trimming and everything downstream.
# Old outputs are moved aside (not deleted) so earlier results are recoverable.
if [ "${RERUN_FROM_TRIM:-0}" = "1" ]; then
	ARCHIVE="${WD}/${CURRENT_HOG}/old_trim_$(date +%y%m%d_%H%M%S)"
	mkdir -p "$ARCHIVE"
	for f in \
		"$CLIPKIT_DIR" \
		"${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.fasta" \
		"${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.dedup.nex" \
		"${WD}/${CURRENT_HOG}/iqtree" \
		"${WD}/${CURRENT_HOG}/${CURRENT_HOG}.orb_fg.tree" \
		"${WD}/${CURRENT_HOG}/${CURRENT_HOG}.nonorb_fg.tree"; do
		[ -e "$f" ] && mv "$f" "$ARCHIVE/"
	done
	# remove-duplicates may also write sidecar files next to the dedup nexus
	mv "${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.dedup."* "$ARCHIVE/" 2>/dev/null || true
	echo "Archived previous trimming-and-later outputs to ${ARCHIVE}"
fi

# Runs clipkit and records its stdout. Args: <mode label> <output prefix> <clipkit mode args...>
run_clipkit() {
	local label="$1" prefix="$2"
	shift 2
	local out="${CLIPKIT_DIR}/${prefix}.fasta"
	local stdout_file="${CLIPKIT_DIR}/${prefix}.stdout.txt"
	clipkit "$MACSE_NT_FILE" "$@" \
		--codon --sequence_type nt --remove_stop_codons all \
		--log --plot_trim_report \
		--output "$out" > "$stdout_file"
	cat "$stdout_file"
	CK_OUT="$out"
	CK_GAPS="$(awk -F': ' '/^Gaps threshold:/ {print $2}' "$stdout_file")"
	# Append stats row (one file per HOG so parallel array tasks don't collide)
	awk -F': ' -v hog="$CURRENT_HOG" -v mode="$label" -v gaps="$CK_GAPS" '
		/^Original length:/ {orig=$2}
		/^Number of sites kept:/ {kept=$2}
		/^Number of sites trimmed:/ {trimmed=$2}
		/^Percentage of alignment trimmed:/ {pct=$2; sub(/%/, "", pct)}
		END {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", hog, mode, gaps, orig, kept, trimmed, pct}
	' "$stdout_file" >> "$CLIPKIT_STATS"
}

if [ -f "$TRIM_FILE" ]; then
	echo "Trimmed file ${TRIM_FILE} exists; on to remove-duplicates."
else
	mkdir -p "$CLIPKIT_DIR"
	printf "hog\tmode\tgaps_threshold\toriginal_length\tsites_kept\tsites_trimmed\tpct_trimmed\n" > "$CLIPKIT_STATS"

	run_clipkit smart-gap "${CURRENT_HOG}_NT.smart-gap" -m smart-gap
	FINAL_CK="$CK_OUT"

	if [ -z "$CK_GAPS" ]; then
		echo "Could not parse ClipKIT gaps threshold from output" >&2
		exit 1
	fi

	if awk -v g="$CK_GAPS" 'BEGIN {exit !(g > 0.9)}'; then
		echo "smart-gap threshold ${CK_GAPS} > 0.9; rerunning with gappy -g 0.9."
		run_clipkit gappy90 "${CURRENT_HOG}_NT.gappy90" -m gappy -g 0.9
		FINAL_CK="$CK_OUT"
	fi

	# Stable name for downstream steps
	cp "$FINAL_CK" "$TRIM_FILE"
	echo "ClipKIT completed; using ${FINAL_CK} as ${TRIM_FILE}"
fi

# Set STOP_AFTER_TRIM=1 to stop once trimming is done (skips dedup, IQ-TREE, labelling)
if [ "${STOP_AFTER_TRIM:-0}" = "1" ]; then
	echo "STOP_AFTER_TRIM=1; stopping after trimming."
	exit 0
fi

#######################
## REMOVE-DUPLICATES ##
#######################

DEDUP_FILE="${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.dedup.nex"

# Check if remove-duplicates has already completed for this HOG
if [ -f "$DEDUP_FILE" ]; then
	echo "Deduplicated file ${DEDUP_FILE} exists; on to IQ-TREE."
else
	# Run remove-duplicates
	hyphy "${HYPHY_ANALYSES}/remove-duplicates/remove-duplicates.bf" \
		"CPU=${SLURM_CPUS_PER_TASK}" \
		--msa "$TRIM_FILE" \
		--output "$DEDUP_FILE" \
		ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

############
## IQTREE ##
############

IQTREE_FILE="${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}.treefile"

#IQ-TREE errors on a finished checkpoint without -redo; a treefile alone doesn't prove the run finished
mkdir -p "${WD}/${CURRENT_HOG}/iqtree/"
IQTREE_LOG="${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}.log"
if [ -f "$IQTREE_LOG" ] && grep -q "Total wall-clock time used" "$IQTREE_LOG"; then
	echo "IQ-TREE run for ${CURRENT_HOG} already finished; on to label-tree."
elif ! iqtree -s "$DEDUP_FILE" \
	-m MFP \
	-T "$SLURM_CPUS_PER_TASK" \
	--prefix "${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}"; then
	# FreeRate EM can abort (ratefree.cpp assertion); retry without +R models
	echo "IQ-TREE failed; retrying without FreeRate (R) models."
	iqtree -s "$DEDUP_FILE" \
		-m MFP \
		-mrate E,I,G \
		-T "$SLURM_CPUS_PER_TASK" \
		-redo \
		--prefix "${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}"
fi

################
## LABEL-TREE ##
################

ORB_LIST="${REPO_ROOT}/data/orbweavers-list.txt"
NONORB_LIST="${REPO_ROOT}/data/non-orbweavers-list.txt"

#label orb-weavers
ORB_TREE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}.orb_fg.tree"
if [ -f "$ORB_TREE" ]; then
	echo "The tree ${ORB_TREE} already exists."
else
	SPECIES_REGEX="$(awk 'NF { sub(/\r$/, ""); printf "%s%s", sep, $0; sep = "|" }' "$ORB_LIST")"
	REGEX="^(${SPECIES_REGEX})\\|"

	hyphy "${HYPHY_ANALYSES}/LabelTrees/label-tree.bf" \
		--tree "$IQTREE_FILE" \
		--regexp "$REGEX" \
		--output "$ORB_TREE" \
		--internal-nodes "All descendants"

fi

#label non-orbweavers
NONORB_TREE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}.nonorb_fg.tree"
if [ -f "$NONORB_TREE" ]; then
	echo "The tree ${NONORB_TREE} already exists."
else
	SPECIES_REGEX="$(awk 'NF { sub(/\r$/, ""); printf "%s%s", sep, $0; sep = "|" }' "$NONORB_LIST")"
	REGEX="^(${SPECIES_REGEX})\\|"

	hyphy "${HYPHY_ANALYSES}/LabelTrees/label-tree.bf" \
		--tree "$IQTREE_FILE" \
		--regexp "$REGEX" \
		--output "$NONORB_TREE" \
		--internal-nodes "All descendants"

fi
