#!/bin/bash

#SBATCH --job-name=260304_get_nuc_seqs
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=00:05:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-3375
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

#run the program from desired CDS output directory

module load anaconda
conda activate /data/agordus1/conda_envs/seqkit/

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

HOG_LIST=/home/crunnel2/kelvin-scratch/assets/N5_occ48_remaining.txt
NX_TSV=/home/crunnel2/orb-selection/assets/N5.tsv
TD_DIR=/data/agordus1/crunnel2/transdecoder

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

python /home/crunnel2/kelvin-scratch/scripts/rockfish/get_nuc_seqs_slurm.py $NX_TSV $CURRENT_HOG $TD_DIR

conda deactivate