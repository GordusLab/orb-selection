#!/bin/bash

# Purpose: run TransDecoder.Predict using retained blastp hits.

# Usage: scripts/01_pre_processing/process_fsas.TD-P.sh <input_dir> <cd_hit_fasta_dir> <blastp_dir> <td_p_complete_dir>

INPUT_DIR=$1
CD_HIT_FASTA_DIR=$2
BLASTP_DIR=$3
TD_P_COMPLETE_DIR=$4

mkdir -p "$TD_P_COMPLETE_DIR"

for dir in "$INPUT_DIR"/*; do
	[ -d "$dir" ] || continue
	NAME=$(basename "$dir")
	cd "$dir" || exit 1
	time TransDecoder.Predict -t "$CD_HIT_FASTA_DIR/$NAME.cd-hit-est" \
		--retain_blastp_hits "$BLASTP_DIR/$NAME.cd-hit-est.blastp.out"
	mkdir -p "$TD_P_COMPLETE_DIR/$NAME"
	mv "$NAME".* "$TD_P_COMPLETE_DIR/$NAME/" 2>/dev/null
	mv pipeliner.* "$TD_P_COMPLETE_DIR/$NAME/" 2>/dev/null
done
