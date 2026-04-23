#!/bin/bash
#SBATCH --job-name=annotate_hogs
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-10000%100 # Adjust array size to number of HOGs and limit parallel jobs

# Define paths
HOG_LIST_FILE="data/N5_hogs_min_occ_30.txt"
N0_TSV="data/N5.tsv"
ORTHO_DIR="data/orthofinder/Results_Apr16/Orthogroup_Sequences"
REF_DB="data/blast_db/dmel_db"
OUTPUT_DIR="results/annotated_hogs"
SCRIPT_PATH="scripts/05_enrichment/annotate_one_hog.py"

# Get the HOG for the current array task
HOG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$HOG_LIST_FILE")

# Run the annotation script for the single HOG
echo "Processing HOG: $HOG"
python3 "$SCRIPT_PATH" "$HOG" "$N0_TSV" "$ORTHO_DIR" "$REF_DB" "$OUTPUT_DIR"
