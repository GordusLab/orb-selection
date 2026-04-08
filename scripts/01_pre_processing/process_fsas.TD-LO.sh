#!/bin/bash

# Purpose: run TransDecoder.LongOrfs on each cd-hit transcript FASTA.

INPUT_DIR=${INPUT_DIR:-PATH_TO_CD_HIT_INPUT_DIR}
SOURCE_FASTA_DIR=${SOURCE_FASTA_DIR:-PATH_TO_CD_HIT_FASTA_DIR}
TD_LO_DIR=${TD_LO_DIR:-PATH_TO_TD_LO_OUTPUT_DIR}

mkdir -p "$TD_LO_DIR"

for file in "$INPUT_DIR"/*.cd-hit-est; do
	[ -e "$file" ] || continue
	NAME=$(basename "$file" .cd-hit-est)
	mkdir -p "$TD_LO_DIR/$NAME"
	cd "$TD_LO_DIR/$NAME" || exit 1
	TransDecoder.LongOrfs -t "$SOURCE_FASTA_DIR/$NAME.cd-hit-est"
done
