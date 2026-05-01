#!/bin/bash
#SBATCH --job-name=260423_annotate_hogs
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --array=1-10 
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out


module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection
module load Biopython/1.81-foss-2022b
module load py-numpy/1.24.4

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

# Define paths
HOG_LIST_FILE="/home/crunnel2/orb-selection/data/N5_hogs_min_occ_30.txt"
N0_TSV="/home/crunnel2/orb-selection/data/N5.tsv"
ORTHO_DIR="/scratch4/agordus1/crunnel2/orthofinder_fastas"
REF_DB="/scratch4/agordus1/crunnel2/GCF_043381705.1/cds_from_genomic.fna"
OUTPUT_DIR="/scratch4/agordus1/crunnel2/annotated_hogs"
SCRIPT_PATH="/home/crunnel2/orb-selection/scripts/05_enrichment/annotate_one_hog.py"

# Get the HOG for the current array task
HOG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$HOG_LIST_FILE")

# Run the annotation script for the single HOG
echo "Processing HOG: $HOG"
python3 "$SCRIPT_PATH" "$HOG" "$N0_TSV" "$ORTHO_DIR" "$REF_DB" "$OUTPUT_DIR"
