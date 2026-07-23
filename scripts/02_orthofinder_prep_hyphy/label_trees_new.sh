#!/bin/bash

#SBATCH --job-name=260722_label_trees
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=01:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-4756
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

WD=/data/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt
FG_LIST=/home/crunnel2/orb-selection/data/non-orb-weavers-list.txt
FG_NAME=non_orb_fg

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)

FLTRD_TREE=${WD}/${CURRENT_HOG}.fltrd.tree

#make a copy of the filtered tree to iteratively label
LBLD_TREE=${WD}/${CURRENT_HOG}.${FG_NAME}.tree
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

	#Add a semicolon to the end of the last line of the tree file to avoid errors in HyPhy
	sed -i '$s/$/;/' ${LBLD_TREE}

fi

conda deactivate
