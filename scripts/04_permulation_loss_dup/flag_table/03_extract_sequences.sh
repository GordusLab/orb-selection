#!/bin/bash
# Step 3 of the flag table: extract the proteins and coding sequences of the candidate genes and their partners.
#
# usage: 03_extract_sequences.sh CANDIDATE_DIR ALL_PEP_FAA CDS_DIR OUT_DIR
#   CANDIDATE_DIR  output of 02_select_candidates.py (protein_ids.txt and cds_ids/<species>.txt)
#   ALL_PEP_FAA    all transcriptome proteins, FASTA IDs '<species>__<gene>' (PATH_TO_ALL_PEP_FAA)
#   CDS_DIR        stage 01 TransDecoder output, one '<species>.cd-hit-est.transdecoder.cds' per species (PATH_TO_TRANSDECODER_CDS)
#   OUT_DIR        writes flag_proteins.faa and flag_cds.fa (CDS IDs renamed to '<species>__<gene>')
set -euo pipefail
CAND=$1; ALL_PEP=$2; CDS_DIR=$3; OUT=$4
mkdir -p "$OUT"

seqkit grep -j 1 -f "$CAND/protein_ids.txt" "$ALL_PEP" > "$OUT/flag_proteins.faa"
echo "$(grep -c '>' "$OUT/flag_proteins.faa") proteins"

: > "$OUT/flag_cds.fa"
for f in "$CAND"/cds_ids/*.txt; do
  s=$(basename "$f" .txt)
  seqkit grep -j 1 -f "$f" "$CDS_DIR/$s.cd-hit-est.transdecoder.cds" \
    | seqkit replace -j 1 -p "^(\S+).*" -r "${s}__\$1" >> "$OUT/flag_cds.fa"
done
echo "$(grep -c '>' "$OUT/flag_cds.fa") coding sequences"
