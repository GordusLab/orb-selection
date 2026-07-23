#!/bin/bash

# Runs RELAX for one HOG per SLURM array task.

#SBATCH --job-name=260723_relax
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-4756
#SBATCH -n 3
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

#check if RELAX has already completed for this HOG 
RELAX_OUT=${WD}/${CURRENT_HOG}_RELAX.json

if grep -q "p-value" ${RELAX_OUT}; then
	echo "RELAX already complete."
else
	#run relax
	hyphy relax \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}.fltrd.fasta \
	 --tree ${WD}/${CURRENT_HOG}.orb_fg.tree \
	 --test "Unlabeled branches" \
	 --multiple-hits Double+Triple \
	 --models Minimal \
	 --srv Yes \
 	 ENV="TOLERATE_NUMERICAL_ERRORS=1;" \
	 --output ${RELAX_OUT}
fi

conda deactivate
