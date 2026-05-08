#!/bin/bash
#SBATCH --job-name=260506_phyloglm_dup_obs
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

# module load GCC/11.3.0 GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 foss/2022a foss/2022b foss/2023a gfbf/2022b gfbf/2023a gfbf/2023b intel/2019a intel/2020a intel/2023a intel/2023b intel/2024a
# module load R/4.5.1-gfbf-2023b
export R_LIBS_USER=/home/crunnel2/R/x86_64-pc-linux-gnu-library/4.5

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/

# Set the number of cores for R
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript /home/crunnel2/orb-selection/scripts/04_permulation_loss_dup/olm_dup.R