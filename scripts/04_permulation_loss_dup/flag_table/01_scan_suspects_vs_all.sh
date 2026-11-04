#!/bin/bash
# Step 1 of the flag table: DIAMOND scan of the suspect samples' proteins against all transcriptome proteins.
#
# usage: 01_scan_suspects_vs_all.sh SUSPECTS_FILE ALL_PEP_FAA OUT_DIR [THREADS]
#   SUSPECTS_FILE  suspect species, one per line (for example data/hub_species_to_exclude.txt)
#   ALL_PEP_FAA    all transcriptome proteins, FASTA IDs '<species>__<gene>' (102 species; PATH_TO_ALL_PEP_FAA)
#   OUT_DIR        output folder (PATH_TO_FLAG_WORKDIR)
#
# Pass 1 runs with DIAMOND's default low-complexity masking; suspect proteins with no hit at all (masking hides
# Gly/Ala-rich proteins) are re-scanned with masking off. Both tables go to 02_select_candidates.py.
# Run on a server with nice; on a shared machine keep THREADS small.
set -euo pipefail
SUSPECTS=$1; ALL_PEP=$2; OUT=$3; THREADS=${4:-3}
mkdir -p "$OUT"

[ -f "$OUT/allpep.dmnd" ] || diamond makedb --in "$ALL_PEP" -d "$OUT/allpep"

PATTERN=$(paste -sd"|" "$SUSPECTS")
seqkit grep -j 1 -r -p "^($PATTERN)__" "$ALL_PEP" > "$OUT/suspects.faa"
echo "$(grep -c '>' "$OUT/suspects.faa") suspect proteins"

FMT="qseqid sseqid pident length qcovhsp scovhsp evalue bitscore"
nice -n 19 diamond blastp --threads "$THREADS" --faster -q "$OUT/suspects.faa" -d "$OUT/allpep" \
  --id 97 --query-cover 70 -e 1e-10 -k 40 --max-hsps 1 -f 6 $FMT -o "$OUT/suspects_vs_all.tsv"

# queries with no row at all in pass 1
seqkit seq -n -i "$OUT/suspects.faa" | sort > "$OUT/all_query_ids.txt"
cut -f1 "$OUT/suspects_vs_all.tsv" | sort -u > "$OUT/queries_with_hits.txt"
comm -23 "$OUT/all_query_ids.txt" "$OUT/queries_with_hits.txt" > "$OUT/unaligned_ids.txt"
seqkit grep -j 1 -f "$OUT/unaligned_ids.txt" "$OUT/suspects.faa" > "$OUT/unaligned.faa"
echo "$(wc -l < "$OUT/unaligned_ids.txt") queries had no hit; re-scanning with masking off"

nice -n 19 diamond blastp --threads "$THREADS" --faster --masking 0 -q "$OUT/unaligned.faa" -d "$OUT/allpep" \
  --id 97 --query-cover 70 -e 1e-10 -k 40 --max-hsps 1 -f 6 $FMT -o "$OUT/unaligned_vs_all.tsv"
