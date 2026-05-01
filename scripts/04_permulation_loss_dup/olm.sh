#!/bin/bash
#SBATCH --job-name=260501_phyloglm
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --mailtype=ALL
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/

module load R

# Set the number of cores for R
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript /path/to/orb-selection/scripts/04_permulation_loss_dup/olm.R