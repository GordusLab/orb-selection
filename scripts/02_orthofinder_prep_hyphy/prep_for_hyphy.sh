#!/bin/bash

# Prepares HOG alignments/trees for HyPhy analyses (PREQUAL -> MACSE -> IQ-TREE -> BUSTED -> error-filter).

#SBATCH --array=1-4756
#SBATCH -n 12
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

#directory containing fasta files of cds to be analyzed
HOG_CDS_DIR=PATH_TO_HOG_CDS_DIR

#list of HOG IDs from directory of fasta files
HOG_LIST=${HOG_LIST:-$REPO_ROOT/data/N5.udiv.o75_list.txt}
MACSE_JAR=PATH_TO_MACSE_JARFILE

CURRENT_HOG=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $HOG_LIST)
mkdir -p ${WD}/${CURRENT_HOG}/

#cds file needs to be in prequal folder so prequal will output to that folder
CDS_FILE=${WD}/${CURRENT_HOG}/prequal/${CURRENT_HOG}.fasta

#check if the cds file has already been moved to the prequal directory
if [ -f ${CDS_FILE} ]; then
	echo "CDS file ${CDS_FILE} already moved; on to prequal."
else
	mkdir ${WD}/${CURRENT_HOG}/prequal/
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
MACSE_FILE=${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_NT.fasta

if [ -f ${MACSE_FILE} ]; then
	echo "Macse file ${MACSE_FILE} exists; on to iqtree."
else
	#run macse
	mkdir ${WD}/${CURRENT_HOG}/macse/
	java -jar -Xmx47G "$MACSE_JAR" \
	 -prog alignSequences \
	 -seq ${PREQUAL_FILE} \
	 -out_NT ${MACSE_FILE} \
	 -out_AA ${WD}/${CURRENT_HOG}/macse/${CURRENT_HOG}_AA.fasta
fi

############
## IQTREE ##
############

#check if iqtree has already completed for this HOG
IQTREE_FILE=${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}.treefile

if [ -f ${IQTREE_FILE} ]; then
	echo "IQtree file ${IQTREE_FILE} exists; on to BUSTED."
else
	#run iqtree
	mkdir ${WD}/${CURRENT_HOG}/iqtree/
	iqtree -s ${MACSE_FILE} \
	-af fasta \
	-T AUTO \
	--prefix ${WD}/${CURRENT_HOG}/iqtree/${CURRENT_HOG}
fi

############
## BUSTED ##
############

#run busted without a foreground in order to extract the error-filtered alignment

#could do remove-duplicates here... but I don't think I want to

#check if BUSTED has already completed for this HOG 
BUSTED_LOG=${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_BUSTED.log

if grep -q "**p =" ${BUSTED_LOG}; then
	echo "BUSTED complete; on to error-filter."
else
	#run busted
	mkdir ${WD}/${CURRENT_HOG}/busted/
	time hyphy busted \
	 CPU=12 \
	 --alignment ${MACSE_FILE} \
	 --tree ${IQTREE_FILE} \
	 --multiple-hits Double+Triple \
	 --output ${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_BUSTED.json \
	 --error-sink Yes \
	 > ${BUSTED_LOG}
fi

##################
## ERROR-FILTER ##
##################

#check if error-filter has already completed for this HOG
FLTRD_FASTA=${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.fasta

if [ -f ${FLTRD_FASTA} ]; then
	echo "Error-filtered file ${FLTRD_FASTA} exists; HOG is ready for HYPHY."
else
	EF_OUT=${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_error-fltrd.nxh

	#extract error-filtered alignment from busted results
	hyphy error-filter \
	 CPU=12 \
	 --json ${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_BUSTED.json \
	 --output ${EF_OUT} \
	 --output-json ${WD}/${CURRENT_HOG}/busted/${CURRENT_HOG}_error-fltrd.json

	#split up tree and alignment from error-filter results
	if grep -q "(" ${EF_OUT}; then
		grep "(" ${EF_OUT} > ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.tree
		sed -i '$s/$/;/' ${WD}/${CURRENT_HOG}/${CURRENT_HOG}.fltrd.tree
		grep -v "(" ${EF_OUT} > ${FLTRD_FASTA}
	else
		echo "Error-filtering failed; check results."
	fi
fi

conda deactivate