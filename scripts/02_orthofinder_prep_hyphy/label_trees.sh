#!/bin/bash

# Labels foreground branches on filtered trees using HyPhy LabelTrees.
# Arguments: <CONDA_ENV_PATH> <REPORT_DIR> <WORK_DIR> <HYPHY_ANALYSES_DIR> <FG_ORB_OR_NONORB>

#SBATCH --array=1-4756
#SBATCH --output=reports/%x/%A_%a.out
#SBATCH --error=reports/%x/%A_%a.err

module load anaconda
CONDA_ENV_PATH=$1
conda activate "$CONDA_ENV_PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

#make directory to store slurm reports
REPORT_DIR=$2
mkdir -p "$REPORT_DIR/$SBATCH_JOB_NAME/"

WD=$3
HOG_LIST=${HOG_LIST:-$REPO_ROOT/data/N5.udiv.o75_list.txt}

if [ "$5" == "orb" ]; then
	FG_LIST=${FG_LIST:-$REPO_ROOT/data/orbweavers-list.txt}
	FG_NAME="orb_fg"
elif [ "$5" == "nonorb" ]; then
	FG_LIST=${FG_LIST:-$REPO_ROOT/data/nonorbweavers-list.txt}
	FG_NAME="nonorb_fg"
else
	echo "Error: Invalid foreground group specified. Use 'orb' or 'nonorb'."
	exit 1
fi

#path to hyphy analyses repo clone: https://github.com/veg/hyphy-analyses/
HYPHY_ANALYSES_DIR=$4

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
		hyphy "$HYPHY_ANALYSES_DIR/LabelTrees/label-tree.bf" \
		 --tree ${LBLD_TREE} \
		 --regexp $REGEX \
		 --output ${LBLD_TREE} 
	done < ${FG_LIST}

	sed -i '$s/$/;/' ${LBLD_TREE}
fi

conda deactivate
