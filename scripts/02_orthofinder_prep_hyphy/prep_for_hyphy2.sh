#!/bin/bash

# Prepares HOG alignments/trees for HyPhy analyses.
# Usage: scripts/02_orthofinder_prep_hyphy/prep_for_hyphy2.sh <conda_env_path> <work_dir> <hog_cds_dir> <macse_jar> <hyphy_analyses>
# Note: The script assumes that the input HOG CDS files are in FASTA format and that MACSE is available at the specified JAR path.

# RUN FROM REPO ROOT
# sbatch scripts/02_orthofinder_prep_hyphy/prep_for_hyphy2.sh ~/anaconda3/envs/hyphy-new/ ~/scratch/hyphy_wd_260929 ~/scratch/hyphy_wd_260929/HOG_CDS/ ~/bin/macse_v2.07.jar ~/bin/hyphy-analyses/

#SBATCH --job-name=261001_prep_for_hyphy2
#SBATCH --partition=shared
#SBATCH --account=agordus1
#SBATCH --time=12:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out
#SBATCH --error=/data/agordus1/crunnel2/reports/%x/%A_%a.err

set -e

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

module load anaconda
CONDA_ENV_PATH=$1
conda activate "$CONDA_ENV_PATH"

REPO_ROOT="$(pwd)"

WD=$2

#directory containing fasta files of cds to be analyzed
HOG_CDS_DIR=$3

#list of HOG IDs from directory of fasta files
HOG_LIST=${REPO_ROOT}/data/N5.udiv.o75_list.txt
MACSE_JAR=$4

HYPHY_ANALYSES=$5

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)
mkdir -p ${WD}/${CURRENT_HOG}/

#cds file needs to be in prequal folder so prequal will output to that folder
CDS_FILE=${WD}/${CURRENT_HOG}/prequal/${CURRENT_HOG}.fasta

#check if the cds file has already been moved to the prequal directory
if [ -f ${CDS_FILE} ]; then
	echo "CDS file ${CDS_FILE} already moved; on to prequal."
else
	mkdir -p ${WD}/${CURRENT_HOG}/prequal/
	mv ${HOG_CDS_DIR}/${CURRENT_HOG}.fasta ${WD}/${CURRENT_HOG}/prequal/
fi

#############
## PREQUAL ##
#############

#check if prequal has already completed for this HOG
PREQUAL_FILE=${WD}/${CURRENT_HOG}/prequal/${CURRENT_HOG}.fasta.dna.filtered

if [ -f ${PREQUAL_FILE} ]; then
	echo "Prequal file ${PREQUAL_FILE} exists; on to macse."
else
	#run prequal
	prequal -dosummary ${CDS_FILE}
fi

###########
## MACSE ##
###########

#check if macse has already completed for this HOG
MACSE_FILE=${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_AA.fasta

if [ -f ${MACSE_FILE} ]; then
	echo "Macse file ${MACSE_FILE} exists; on to TrimAl."
else
	#run macse
	mkdir -p ${WD}/${CURRENT_HOG}/macse/
	java -jar -Xmx47G "$MACSE_JAR" \
	 -prog alignSequences \
	 -seq ${PREQUAL_FILE} \
	 -out_NT ${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.fasta \
	 -out_AA ${MACSE_FILE}
fi

############
## TRIMAL ##
############

# Trim protein alignments using trimAl, then back-translate to nucleotide sequences
TRIMAL_FILE=${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.fasta

PREQUAL_CDS=${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_CDS_for_trimAl.fasta
# replace all X/x with N/n in the prequal file for trimAl back-translation
sed '/^>/! s/X/N/g; /^>/! s/x/n/g' "$PREQUAL_FILE" > "$PREQUAL_CDS"

# check if trimAl has already completed for this HOG
if [ -f ${TRIMAL_FILE} ]; then
	echo "TrimAl file ${TRIMAL_FILE} exists; on to remove-duplicates."
else
	#run trimAl
	trimal -in ${MACSE_FILE} -backtrans "$PREQUAL_CDS" -out ${TRIMAL_FILE} -gappyout
	echo "TrimAl completed for ${TRIMAL_FILE}"
fi

#######################
## REMOVE-DUPLICATES ##
#######################

DEDUP_FILE=${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.trim.dedup.fasta

# check if remove-duplicates has already completed for this HOG
if [ -f ${DEDUP_FILE} ]; then
	echo "Deduplicated file ${DEDUP_FILE} exists; on to IQ-TREE."
else
	#run remove-duplicates
	hyphy ${HYPHY_ANALYSES}/remove-duplicates/remove-duplicates.bf \
	 CPU=${SLURM_NTASKS} \
	 --msa ${TRIMAL_FILE} \
	 --output ${DEDUP_FILE} \
	 ENV="TOLERATE_NUMERICAL_ERRORS=1;"
fi

############
## IQTREE ##
############

#check if iqtree has already completed for this HOG
IQTREE_FILE=${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}.treefile

if [ -f ${IQTREE_FILE} ]; then
	echo "IQtree file ${IQTREE_FILE} exists; on to LabelTrees."
else
	#run iqtree
	mkdir -p ${WD}/${CURRENT_HOG}/iqtree/
	iqtree -s ${DEDUP_FILE} \
	-m MFP \
	-ntmax 3 \
	-T AUTO \
	--prefix ${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}
fi

################
## LABEL-TREE ##
################

ORB_LIST=${REPO_ROOT}/data/orbweavers-list.txt
NONORB_LIST=${REPO_ROOT}/data/non-orbweavers-list.txt

#label orb-weavers
ORB_TREE=${WD}/${CURRENT_HOG}/${CURRENT_HOG}.orb_fg.parsimony.tree
if [ -f "${ORB_TREE}" ]; then
	echo "The tree ${ORB_TREE} already exists."
else
	SPECIES_REGEX=$(awk 'NF { sub(/\r$/, ""); printf "%s%s", sep, $0; sep = "|" }' "$ORB_LIST")
	REGEX="^(${SPECIES_REGEX})\\|"

	hyphy "$HYPHY_ANALYSES/LabelTrees/label-tree.bf" \
	--tree "$IQTREE_FILE" \
	--regexp "$REGEX" \
	--output "$ORB_TREE" \
	--internal-nodes "Parsimony"

fi

#label non-orbweavers
NONORB_TREE=${WD}/${CURRENT_HOG}/${CURRENT_HOG}.nonorb_fg.parsimony.tree
if [ -f "${NONORB_TREE}" ]; then
	echo "The tree ${NONORB_TREE} already exists."
else
	SPECIES_REGEX=$(awk 'NF { sub(/\r$/, ""); printf "%s%s", sep, $0; sep = "|" }' "$NONORB_LIST")
	REGEX="^(${SPECIES_REGEX})\\|"

	hyphy "$HYPHY_ANALYSES/LabelTrees/label-tree.bf" \
	--tree "$IQTREE_FILE" \
	--regexp "$REGEX" \
	--output "$NONORB_TREE" \
	--internal-nodes "Parsimony"

fi

# ############
# ## BUSTED ##
# ############

# #run busted without a foreground in order to extract the error-filtered alignment

# #could do remove-duplicates here... but I don't think I want to

# #check if BUSTED has already completed for this HOG 
# BUSTED_LOG=${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_BUSTED.log

# if grep -q "**p =" ${BUSTED_LOG}; then
# 	echo "BUSTED complete; on to error-filter."
# else
# 	#run busted
# 	mkdir -p ${WD}/${CURRENT_HOG}/busted/
# 	time hyphy busted \
# 	 CPU=12 \
# 	 --alignment ${MACSE_FILE} \
# 	 --tree ${IQTREE_FILE} \
# 	 --multiple-hits Double+Triple \
# 	 --output ${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_BUSTED.json \
# 	 --error-sink Yes \
# 	 > ${BUSTED_LOG}
# fi

# ##################
# ## ERROR-FILTER ##
# ##################

# #check if error-filter has already completed for this HOG
# FLTRD_FASTA=${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.fasta

# if [ -f ${FLTRD_FASTA} ]; then
# 	echo "Error-filtered file ${FLTRD_FASTA} exists; HOG is ready for HYPHY."
# else
# 	EF_OUT=${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_error-fltrd.nxh

# 	#extract error-filtered alignment from busted results
# 	hyphy error-filter \
# 	 CPU=12 \
# 	 --json ${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_BUSTED.json \
# 	 --output ${EF_OUT} \
# 	 --output-json ${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_error-fltrd.json

# 	#split up tree and alignment from error-filter results
# 	if grep -q "(" ${EF_OUT}; then
# 		grep "(" ${EF_OUT} > ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.tree
# 		sed -i '$s/$/;/' ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.tree
# 		grep -v "(" ${EF_OUT} > ${FLTRD_FASTA}
# 	else
# 		echo "Error-filtering failed; check results."
# 	fi
# fi

# conda deactivate