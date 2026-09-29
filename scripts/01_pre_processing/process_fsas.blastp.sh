#!/bin/bash

# Purpose: run blastp on longest TransDecoder ORFs and save tabular hits.
# Usage: scripts/01_pre_processing/process_fsas.blastp.sh <input_dir> <td_lo_to_blastp_dir> <blastp_out_dir> <uniprot_db> <td_p_run_dir>

INPUT_DIR=$1
TD_LO_TO_BLASTP_DIR=$2
BLASTP_OUT_DIR=$3
UNIPROT_DB=$4
TD_P_RUN_DIR=$5

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
