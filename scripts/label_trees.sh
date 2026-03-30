#!/bin/bash

#SBATCH --job-name=250605_label_trees
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=00:10:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-4756
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

WD=/data/agordus1/crunnel2/060125_N5_o75
HOG_LIST=${WD}/HOG_CDS/N5.udiv.o75_list.txt
FG_LIST=/home/crunnel2/kelvin-scratch/data/family-tree/non-orbweavers-list.txt
FG_NAME=nonorb_fg

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

FLTRD_TREE=${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.tree

#make a copy of the filtered tree to iteratively label
LBLD_TREE=${WD}/${CURRENT_HOG}/${CURRENT_HOG}.${FG_NAME}.tree
if [ -f ${LBLD_TREE} ]; then
	echo "The tree ${LBLD_TREE} already exists."
else
	cp ${FLTRD_TREE} ${LBLD_TREE}

	while read p; do
		REGEX="^${p}"
		hyphy /home/crunnel2/bin/hyphy-analyses/LabelTrees/label-tree.bf \
		 --tree ${LBLD_TREE} \
		 --regexp $REGEX \
		 --output ${LBLD_TREE} 
	done < ${FG_LIST}

	sed -i '$s/$/;/' ${LBLD_TREE}
fi

conda deactivate
