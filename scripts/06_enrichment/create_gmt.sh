#!/usr/bin/env bash

set -euo pipefail

gaf="${1:-data/udiv_gene_ontology.gaf}"
local_gmt="${2:-data/GO.gmt}"
output="${3:-data/udiv_GO.gmt}"

work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

awk -F '\t' '!/^!/ {print $5 "\t" $2}' "$gaf" | sort -u > "$work_dir/pairs.tsv"
cut -f1 "$work_dir/pairs.tsv" | sort -u > "$work_dir/terms.txt"

awk -F '\t' '!seen[$1]++ {print $1 "\t" $2}' "$local_gmt" | sort -u > "$work_dir/local_descriptions.tsv"
comm -23 "$work_dir/terms.txt" <(cut -f1 "$work_dir/local_descriptions.tsv") > "$work_dir/missing_terms.txt"

: > "$work_dir/online_descriptions.tsv"
while IFS= read -r term; do
    label=$(curl -L --fail --silent --show-error "https://api.geneontology.org/api/ontology/term/$term" | jq -r '.label // empty')
    if [[ -n "$label" ]]; then
        printf '%s\t%s\n' "$term" "$label" >> "$work_dir/online_descriptions.tsv"
    else
        printf '%s\t[description unavailable]\n' "$term" >> "$work_dir/online_descriptions.tsv"
    fi
done < "$work_dir/missing_terms.txt"

cat "$work_dir/local_descriptions.tsv" "$work_dir/online_descriptions.tsv" | sort -u > "$work_dir/descriptions.tsv"

awk -F '\t' '
    NR == FNR {description[$1] = $2; next}
    {genes[$1] = genes[$1] "\t" $2}
    END {
        for (term in genes) {
            label = description[term]
            if (label == "") label = "[description unavailable]"
            print term "\t" label genes[term]
        }
    }
' "$work_dir/descriptions.tsv" "$work_dir/pairs.tsv" | sort -t $'\t' -k1,1 > "$output"

printf 'Wrote %s (%s terms, %s unique term-gene pairs; %s online descriptions requested)\n' \
    "$output" "$(wc -l < "$work_dir/terms.txt")" "$(wc -l < "$work_dir/pairs.tsv")" "$(wc -l < "$work_dir/missing_terms.txt")"