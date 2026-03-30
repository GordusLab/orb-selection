#!/bin/bash

#SBATCH --job-name=250218_busted_alphabet_errors
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-2
#SBATCH -n 12
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/${SBATCH_JOB_NAME}/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

WD=/data/agordus1/crunnel2/060125_N5_o75
HOG_LIST=${WD}/HOG_CDS/alphabetbois.txt
FG_NAME=ow_fg_cons	

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

#check if BUSTED-PH has already completed for this HOG 
BUSTEDPH_OUT=${WD}/${CURRENT_HOG}/${CURRENT_HOG}_BUSTED-PH.json

if grep -q "p-value" ${BUSTEDPH_OUT}; then
	echo "BUSTED-PH already complete."
else
	#run busted-ph
	hyphy /home/crunnel2/bin/hyphy-analyses/BUSTED-PH/BUSTED-PH.bf \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.fasta \
	 --tree ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.tree \
	 --branches Foreground \
	 --output ${BUSTEDPH_OUT} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

conda deactivate