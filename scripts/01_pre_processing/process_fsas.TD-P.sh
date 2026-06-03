#!/bin/bash

# Purpose: run TransDecoder.Predict using retained blastp hits.

INPUT_DIR=PATH_TO_TD_LO_OUTPUT_DIR
CD_HIT_FASTA_DIR=PATH_TO_CD_HIT_FASTA_DIR
BLASTP_DIR=PATH_TO_BLASTP_OUTPUT_DIR
TD_P_COMPLETE_DIR=PATH_TO_TD_P_COMPLETE_DIR

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
