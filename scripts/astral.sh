#!/bin/bash

#SBATCH --job-name=250807_astral
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=48
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/astral

astral-pro3 \
  -t ${SLURM_CPUS_PER_TASK} \
  -o astral_speciestree.nw \
  -i multitree.nw \ 
  -g pruned_speciestree.nw \
  2> astral_speciestree.log