#!/bin/bash

# Labels foreground branches on filtered trees using HyPhy LabelTrees.

#SBATCH --array=1-4756
#SBATCH --output=reports/%x/%A_%a.out

module load anaconda
CONDA_ENV_PATH=PATH_TO_CONDA_ENV
conda activate "$CONDA_ENV_PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

#make directory to store slurm reports
REPORT_DIR=PATH_TO_REPORT_DIR
mkdir -p "$REPORT_DIR/$SBATCH_JOB_NAME/"

WD=PATH_TO_EXTERNAL_WORK_DIR
HOG_LIST=${HOG_LIST:-$REPO_ROOT/data/N5.udiv.o75_list.txt}
FG_LIST=${FG_LIST:-$REPO_ROOT/data/non-orbweavers-list.txt}
FG_NAME=nonorb_fg

#path to hyphy analyses repo clone: https://github.com/veg/hyphy-analyses/
HYPHY_ANALYSES_DIR=PATH_TO_HYPHY_ANALYSES_DIR_CLONE

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
		hyphy "$HYPHY_ANALYSES_DIR/LabelTrees/label-tree.bf" \
		 --tree ${LBLD_TREE} \
		 --regexp $REGEX \
		 --output ${LBLD_TREE} 
	done < ${FG_LIST}

	sed -i '$s/$/;/' ${LBLD_TREE}
fi

conda deactivate
