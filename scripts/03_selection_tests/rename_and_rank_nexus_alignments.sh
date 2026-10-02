#!/usr/bin/env bash

# Rename legacy NEXUS alignments from .fasta to .nex and rank by NTAX * NCHAR.
# Usage: rename_and_rank_nexus_alignments.sh [work_dir] [hog_list] [top_n]

set -euo pipefail

WD="${1:-/scratch4/agordus1/crunnel2/hyphy_wd_260929}"
HOG_LIST="${2:-/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt}"
TOP_N="${3:-10}"

if [[ ! -r "$HOG_LIST" ]]; then
    printf 'Cannot read HOG list: %s\n' "$HOG_LIST" >&2
    exit 1
fi

if [[ ! -d "$WD" ]]; then
    printf 'Work directory does not exist: %s\n' "$WD" >&2
    exit 1
fi

if [[ ! "$TOP_N" =~ ^[1-9][0-9]*$ ]]; then
    printf 'top_n must be a positive integer: %s\n' "$TOP_N" >&2
    exit 1
fi

rank_file="$(mktemp)"
trap 'rm -f "$rank_file"' EXIT

array_index=0
renamed=0
already_named=0
missing=0
invalid=0

while IFS= read -r hog || [[ -n "$hog" ]]; do
    array_index=$((array_index + 1))
    hog="${hog%$'\r'}"
    [[ -z "$hog" ]] && continue

    old_file="${WD}/${hog}/macse/${hog}_NT.trim.dedup.fasta"
    nexus_file="${WD}/${hog}/macse/${hog}_NT.trim.dedup.nex"

    if [[ -e "$old_file" && -e "$nexus_file" ]]; then
        printf 'Both old and new filenames exist for %s; refusing to overwrite either.\n' "$hog" >&2
        exit 1
    elif [[ -f "$old_file" ]]; then
        mv -- "$old_file" "$nexus_file"
        renamed=$((renamed + 1))
    elif [[ -f "$nexus_file" ]]; then
        already_named=$((already_named + 1))
    else
        missing=$((missing + 1))
        continue
    fi

    if ! grep -Eiq '^[[:space:]]*#NEXUS([[:space:]]|$)' "$nexus_file"; then
        printf 'Not recognized as NEXUS; skipped ranking: %s\n' "$nexus_file" >&2
        invalid=$((invalid + 1))
        continue
    fi

    dimensions="$(awk '
        function get_dimension(line, key, matched) {
            if (match(line, key "[[:space:]]*=[[:space:]]*[0-9]+")) {
                matched = substr(line, RSTART, RLENGTH)
                sub(/.*=[[:space:]]*/, "", matched)
                return matched + 0
            }
            return 0
        }
        {
            line = toupper($0)
            if (!ntax) ntax = get_dimension(line, "NTAX")
            if (!nchar) nchar = get_dimension(line, "NCHAR")
        }
        END {
            if (ntax > 0 && nchar > 0) printf "%d\t%d\n", ntax, nchar
        }
    ' "$nexus_file")"

    if [[ -z "$dimensions" ]]; then
        printf 'Could not read NTAX and NCHAR dimensions; skipped ranking: %s\n' "$nexus_file" >&2
        invalid=$((invalid + 1))
        continue
    fi

    IFS=$'\t' read -r ntax nchar <<< "$dimensions"
    workload=$((ntax * nchar))
    printf '%s\t%s\t%s\t%s\t%s\n' "$array_index" "$hog" "$ntax" "$nchar" "$workload" >> "$rank_file"
done < "$HOG_LIST"

printf 'Renamed: %s; already .nex: %s; missing: %s; invalid/unranked: %s\n' \
    "$renamed" "$already_named" "$missing" "$invalid" >&2
printf 'Array index\tHOG\tNTAX\tNCHAR\tNTAX*NCHAR (estimate)\n'
sort -t $'\t' -k5,5n "$rank_file" | sed -n "1,${TOP_N}p"