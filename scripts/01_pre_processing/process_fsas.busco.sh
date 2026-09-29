#!/bin/bash

# Purpose: run BUSCO on TransDecoder peptide outputs for each sample.
# Usage: scripts/01_pre_processing/process_fsas.busco.sh <input_dir> <busco_work_dir> <td_p_busco_run_dir> <td_p_busco_complete_dir>

INPUT_DIR=$1
BUSCO_WORK_DIR=$2
TD_P_BUSCO_RUN_DIR=$3
TD_P_BUSCO_COMPLETE_DIR=$4

mkdir -p "$BUSCO_WORK_DIR"
mkdir -p "$TD_P_BUSCO_COMPLETE_DIR"

for dir in "$INPUT_DIR"/*; do
	[ -d "$dir" ] || continue
	NAME=$(basename "$dir")
	mkdir -p "$BUSCO_WORK_DIR/$NAME"
	cd "$BUSCO_WORK_DIR/$NAME" || exit 1
	busco -i "$TD_P_BUSCO_RUN_DIR/$NAME/$NAME.cd-hit-est.transdecoder.pep" \
		-l arachnida \
		-o "$NAME.cd-hit-est.transdecoder.busco" \
		-m prot \
		-c 60 \
		--datasets_version odb10
	mv "$TD_P_BUSCO_RUN_DIR/$NAME" "$TD_P_BUSCO_COMPLETE_DIR/$NAME"
done
