#!/bin/bash

# Purpose: run blastp on longest TransDecoder ORFs and save tabular hits.

INPUT_DIR=PATH_TO_CD_HIT_INPUT_DIR
TD_LO_TO_BLASTP_DIR=PATH_TO_TD_LO_TO_BLASTP_DIR
BLASTP_OUT_DIR=PATH_TO_BLASTP_OUTPUT_DIR
UNIPROT_DB=PATH_TO_UNIPROT_SPROT_FASTA
TD_P_RUN_DIR=PATH_TO_TD_P_RUN_DIR

mkdir -p "$BLASTP_OUT_DIR"
mkdir -p "$TD_P_RUN_DIR"

for file in "$INPUT_DIR"/*.cd-hit-est; do
  [ -e "$file" ] || continue
  NAME=$(basename "$file" .cd-hit-est)
  echo "$NAME"
  time blastp -query "$TD_LO_TO_BLASTP_DIR/$NAME/$NAME.cd-hit-est.transdecoder_dir/longest_orfs.pep" \
	-db "$UNIPROT_DB" \
	-max_target_seqs 1 \
	-outfmt 6 \
	-evalue 1e-5 \
	-num_threads 36 \
	> "$BLASTP_OUT_DIR/$NAME.cd-hit-est.blastp.out"
  mv "$file" "$TD_P_RUN_DIR/"
done
