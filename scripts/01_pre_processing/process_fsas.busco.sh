#!/bin/bash

# Purpose: run BUSCO on TransDecoder peptide outputs for each sample.

INPUT_DIR=${INPUT_DIR:-PATH_TO_TD_P_OUTPUT_PARENT_DIR}
BUSCO_WORK_DIR=${BUSCO_WORK_DIR:-PATH_TO_BUSCO_WORK_DIR}
TD_P_BUSCO_RUN_DIR=${TD_P_BUSCO_RUN_DIR:-PATH_TO_TD_P_BUSCO_RUN_DIR}
TD_P_BUSCO_COMPLETE_DIR=${TD_P_BUSCO_COMPLETE_DIR:-PATH_TO_TD_P_BUSCO_COMPLETE_DIR}

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
