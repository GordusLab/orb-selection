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

WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)
OW_TREE=${WD}/${CURRENT_HOG}.${orb_fg}.tree
NOW_TREE=${WD}/${CURRENT_HOG}.${non_orb_fg}.tree

mv $OW_TREE $OW_TREE.tmp

gotree brlen clear -i $OW_TREE.tmp -o $OW_TREE && rm $OW_TREE.tmp

mv $NOW_TREE $NOW_TREE.tmp

gotree brlen clear -i $NOW_TREE.tmp -o $NOW_TREE && rm $NOW_TREE.tmp