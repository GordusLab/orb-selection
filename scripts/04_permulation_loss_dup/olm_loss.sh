#!/bin/bash
#SBATCH --job-name=260505_phyloglm_loss
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

# module load R/4.5.1-gfbf-2023b
export R_LIBS_USER=/home/crunnel2/R/x86_64-pc-linux-gnu-library/4.5

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/

# Set the number of cores for R
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript /home/crunnel2/orb-selection/scripts/04_permulation_loss_dup/olm_loss.R