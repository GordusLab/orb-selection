#!/bin/bash

#SBATCH --job-name=260820_remove_dups_trim
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=00:15:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=10-4756
#SBATCH -n 1
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda

WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

TRIM_ALGN=${WD}/${CURRENT_HOG}.trim.fasta
if [ -f ${TRIM_ALGN} ]; then
	echo "The trimmed alignment ${TRIM_ALGN} already exists."
else
	conda activate /home/crunnel2/anaconda3/envs/trimal
	trimal -in ${WD}/${CURRENT_HOG}.fltrd.fasta -out ${TRIM_ALGN} -gappyout
	conda deactivate
fi

conda activate /home/crunnel2/anaconda3/envs/hyphy-new

#for orb_fg
hyphy /home/crunnel2/bin/hyphy-analyses/remove-duplicates/remove-duplicates.bf \
 CPU=${SLURM_NTASKS} \
 --msa ${TRIM_ALGN} \
 --tree ${WD}/${CURRENT_HOG}.orb_fg.tree \
 --output ${WD}/${CURRENT_HOG}.orb_fg.trim.nex \
 ENV="TOLERATE_NUMERICAL_ERRORS=1;"

#for non_orb_fg
hyphy /home/crunnel2/bin/hyphy-analyses/remove-duplicates/remove-duplicates.bf \
 CPU=${SLURM_NTASKS} \
 --msa ${TRIM_ALGN} \
 --tree ${WD}/${CURRENT_HOG}.non_orb_fg.tree \
 --output ${WD}/${CURRENT_HOG}.non_orb_fg.trim.nex \
 ENV="TOLERATE_NUMERICAL_ERRORS=1;"


conda deactivate
