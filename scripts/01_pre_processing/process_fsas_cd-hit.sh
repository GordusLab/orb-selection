#!/bin/bash

# Purpose: cluster nucleotide FASTAs with cd-hit-est.
# Usage: scripts/01_pre_processing/process_fsas_cd-hit.sh <input_dir> <cd_hit_out_dir> <cd_hit_complete_dir>

INPUT_DIR=$1
CD_HIT_OUT_DIR=$2
CD_HIT_COMPLETE_DIR=$3

mkdir -p "$CD_HIT_OUT_DIR"
mkdir -p "$CD_HIT_COMPLETE_DIR"

for file in "$INPUT_DIR"/*.fsa_nt; do
	[ -e "$file" ] || continue
	NAME=$(basename "$file" .fsa_nt)
	cd-hit-est -i "$file" -b 25 -M 0 -T 0 -d 128 -p 1 -g 1 -o "$CD_HIT_OUT_DIR/$NAME.cd-hit-est"
	mv "$file" "$CD_HIT_COMPLETE_DIR/"
done
