#!/bin/bash

# Prepares HOG alignments/trees for HyPhy analyses.
# Usage: scripts/02_orthofinder_prep_hyphy/prep_for_hyphy2.sh <conda_env_path> <work_dir> <hog_cds_dir> <macse_jar> <hyphy_analyses>
# Note: The script assumes that the input HOG CDS files are in FASTA format and that MACSE is available at the specified JAR path.
# Assumes that the repo is located at $HOME/orb-selection

# sbatch ~/orb-selection/scripts/02_orthofinder_prep_hyphy/prep_for_hyphy2.sh ~/anaconda3/envs/hyphy-new/ ~/scratch/hyphy_wd_260929 ~/scratch/hyphy_wd_260929/HOG_CDS/ ~/bin/macse_v2.07.jar ~/bin/hyphy-analyses/

#SBATCH --job-name=261006_prep_for_hyphy2_full_iqtree
#SBATCH --partition=shared
#SBATCH --account=agordus1
#SBATCH --time=00:10:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-4756%500
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --output=/scratch4/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/scratch4/agordus1/crunnel2/reports/%x/%A_%a.err

set -e

# Make directory to store SLURM reports
mkdir -p "/scratch4/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/"

CONDA_ENV_PATH="$1"
WD="$2"
HOG_CDS_DIR="$3"
MACSE_JAR="$4"
HYPHY_ANALYSES="$5"

module load anaconda
conda activate "$CONDA_ENV_PATH"
# Cluster modules can leak PYTHONPATH (e.g. a python3.9 biopython) that shadows the env's packages
unset PYTHONPATH

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

# Remove truncated/frameshifted sequences (and mask remaining "!" frameshift codons as gaps), then trim the MACSE codon alignment with
# ClipKIT gappy -g 0.9 (drops only sites that are >90% gaps).
CLIPKIT_DIR="${WD}/${CURRENT_HOG}/clipkit"
FINAL_CK="${CLIPKIT_DIR}/${CURRENT_HOG}_NT.gappy90.fasta"
CLIPKIT_STATS="${CLIPKIT_DIR}/${CURRENT_HOG}_clipkit_stats.tsv"
FILTER_DIR="${WD}/${CURRENT_HOG}/seqfilter"
FILTERED_FILE="${FILTER_DIR}/${CURRENT_HOG}_NT.filtered.fasta"

# Set RERUN_FROM_TRIM=1 to redo trimming and everything downstream.
# Old outputs are moved aside (not deleted) so earlier results are recoverable.
if [ "${RERUN_FROM_TRIM:-0}" = "1" ]; then
	ARCHIVE="${WD}/${CURRENT_HOG}/old_trim_$(date +%y%m%d_%H%M%S)"
	mkdir -p "$ARCHIVE"
	for f in \
		"$CLIPKIT_DIR" \
		"$FILTER_DIR" \
		"${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.fasta" \
		"${WD}/${CURRENT_HOG}/${CURRENT_HOG}_NT.trim.dedup.nex" \
		"${WD}/${CURRENT_HOG}/iqtree" \
		"${WD}/${CURRENT_HOG}/${CURRENT_HOG}.orb_fg.tree" \
		"${WD}/${CURRENT_HOG}/${CURRENT_HOG}.nonorb_fg.tree"; do
		[ -e "$f" ] && mv "$f" "$ARCHIVE/"
	done
	# HyPhy results (BUSTED-PH, RELAX, and their intermediates) depend on the old alignment/tree
	mv "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED"*.json "${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX"*.json "$ARCHIVE/" 2>/dev/null || true
	# remove-duplicates may also write sidecar files next to the dedup nexus
	mv "${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.dedup."* "$ARCHIVE/" 2>/dev/null || true
	echo "Archived previous trimming-and-later outputs to ${ARCHIVE}"
fi

# The final file is only moved into place once ClipKIT succeeds, so its presence means trimming finished
if [ -s "$FINAL_CK" ]; then
	echo "ClipKIT already done; using ${FINAL_CK}; on to remove-duplicates."
else
	mkdir -p "$CLIPKIT_DIR" "$FILTER_DIR"

	# Remove truncated/frameshifted sequences before trimming (log in seqfilter/)
	python "${REPO_ROOT}/scripts/02_orthofinder_prep_hyphy/filter_msa_sequences.py" \
		"$MACSE_NT_FILE" "$FILTERED_FILE" "${FILTER_DIR}/${CURRENT_HOG}_removed.tsv" \
		"${REPO_ROOT}/data/orbweavers-list.txt"

	N_SEQS_IN="$(grep -c '^>' "$MACSE_NT_FILE")"
	N_SEQS_OUT="$(grep -c '^>' "$FILTERED_FILE")"
	N_SEQS_REMOVED=$((N_SEQS_IN - N_SEQS_OUT))

	ck_stdout="$(clipkit "$FILTERED_FILE" -m gappy -g 0.9 \
		--codon --sequence_type nt --remove_stop_codons all \
		-t "$SLURM_CPUS_PER_TASK" \
		--output "${FINAL_CK}.tmp")"
	echo "$ck_stdout"

	# One stats row per HOG (one file per HOG so parallel array tasks don't collide)
	{
		printf "hog\tmode\tgaps_threshold\toriginal_length\tsites_kept\tsites_trimmed\tpct_trimmed\tn_seqs_in\tn_seqs_removed\tn_seqs_out\n"
		awk -F': ' -v hog="$CURRENT_HOG" -v nin="$N_SEQS_IN" -v nrem="$N_SEQS_REMOVED" -v nout="$N_SEQS_OUT" '
			/^Original length:/ {orig=$2}
			/^Number of sites kept:/ {kept=$2}
			/^Number of sites trimmed:/ {trimmed=$2}
			/^Percentage of alignment trimmed:/ {pct=$2; sub(/%/, "", pct)}
			END {printf "%s\tgappy90\t0.9\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", hog, orig, kept, trimmed, pct, nin, nrem, nout}
		' <<< "$ck_stdout"
	} > "$CLIPKIT_STATS"

	if [ ! -s "${FINAL_CK}.tmp" ]; then
		echo "ClipKIT did not produce a non-empty ${FINAL_CK}.tmp" >&2
		exit 1
	fi
	mv "${FINAL_CK}.tmp" "$FINAL_CK"
	echo "ClipKIT completed; using ${FINAL_CK}"
fi

# Set STOP_AFTER_TRIM=1 to stop once trimming is done (skips dedup, IQ-TREE, labelling)
if [ "${STOP_AFTER_TRIM:-0}" = "1" ]; then
	echo "STOP_AFTER_TRIM=1; stopping after trimming."
	exit 0
fi

#######################
## REMOVE-DUPLICATES ##
#######################

DEDUP_FILE="${WD}/${CURRENT_HOG}/${CURRENT_HOG}_NT.trim.dedup.nex"

# Check if remove-duplicates has already completed for this HOG
if [ -f "$DEDUP_FILE" ]; then
	echo "Deduplicated file ${DEDUP_FILE} exists; on to IQ-TREE."
else
	# Run remove-duplicates
	hyphy "${HYPHY_ANALYSES}/remove-duplicates/remove-duplicates.bf" \
		"CPU=${SLURM_CPUS_PER_TASK}" \
		--msa "$FINAL_CK" \
		--output "$DEDUP_FILE" \
		ENV="TOLERATE_NUMERICAL_ERRORS=1;"
	# HyPhy can exit 0 without writing output
	if [ ! -s "$DEDUP_FILE" ]; then
		echo "remove-duplicates did not produce ${DEDUP_FILE}" >&2
		exit 1
	fi
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
