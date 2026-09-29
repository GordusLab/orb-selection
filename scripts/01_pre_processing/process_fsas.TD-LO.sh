#!/bin/bash

# Purpose: run TransDecoder.LongOrfs on each cd-hit transcript FASTA.
# Usage: scripts/01_pre_processing/process_fsas.TD-LO.sh <input_dir> <source_fasta_dir> <td_lo_output_dir>

INPUT_DIR=$1
SOURCE_FASTA_DIR=$2
TD_LO_DIR=$3

mkdir -p "$TD_LO_DIR"

for file in "$INPUT_DIR"/*.cd-hit-est; do
	[ -e "$file" ] || continue
	NAME=$(basename "$file" .cd-hit-est)
	mkdir -p "$TD_LO_DIR/$NAME"
	cd "$TD_LO_DIR/$NAME" || exit 1
	TransDecoder.LongOrfs -t "$SOURCE_FASTA_DIR/$NAME.cd-hit-est"
done
