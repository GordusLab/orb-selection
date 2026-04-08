#!/bin/bash

# Purpose: cluster nucleotide FASTAs with cd-hit-est.

INPUT_DIR=${INPUT_DIR:-PATH_TO_FSA_INPUT_DIR}
CD_HIT_OUT_DIR=${CD_HIT_OUT_DIR:-PATH_TO_CD_HIT_OUTPUT_DIR}
CD_HIT_COMPLETE_DIR=${CD_HIT_COMPLETE_DIR:-PATH_TO_PROCESSED_INPUT_DIR}

mkdir -p "$CD_HIT_OUT_DIR"
mkdir -p "$CD_HIT_COMPLETE_DIR"

for file in "$INPUT_DIR"/*.fsa_nt; do
	[ -e "$file" ] || continue
	NAME=$(basename "$file" .fsa_nt)
	cd-hit-est -i "$file" -b 25 -M 0 -T 0 -d 128 -p 1 -g 1 -o "$CD_HIT_OUT_DIR/$NAME.cd-hit-est"
	mv "$file" "$CD_HIT_COMPLETE_DIR/"
done
