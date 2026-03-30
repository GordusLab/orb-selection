#!/bin/bash

#SBATCH --job-name=250310_relax_final
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 48
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

WD=/data/agordus1/crunnel2/060125_N5_o75
# HOG_LIST=${WD}/HOG_CDS/relax_final_5.txt
FG_NAME=ow_fg_cons

# CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

CURRENT_HOG=N5.HOG0038418

#check if RELAX has already completed for this HOG 
RELAX_OUT=${WD}/${CURRENT_HOG}/${CURRENT_HOG}_RELAX.json

if grep -q "p-value" ${RELAX_OUT}; then
	echo "RELAX already complete."
else
	#run relax
	hyphy relax \
	 CPU=${SLURM_NTASKS} \
	 --alignment ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.fasta \
	 --tree ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.tree \
	 --test "Unlabeled branches" \
	 --multiple-hits Double+Triple \
	 --models Minimal \
	 --srv Yes \
 	 ENV="TOLERATE_NUMERICAL_ERRORS=1;" \
	 --output ${RELAX_OUT}
fi

conda deactivate
