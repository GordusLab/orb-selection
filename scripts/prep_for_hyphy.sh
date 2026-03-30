#!/bin/bash

#SBATCH --job-name=140125_prep_hyphy_max_resources
#SBATCH --partition=parallel
#SBATCH --account=agordus1
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=crunnel2@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-30
#SBATCH -n 12
#SBATCH --output=/data/agordus1/crunnel2/reports/%x/%A_%a.out

module load anaconda
conda activate /home/crunnel2/anaconda3/envs/selection

#make directory to store slurm reports
mkdir -p /data/agordus1/crunnel2/reports/$SBATCH_JOB_NAME/

WD=/data/agordus1/crunnel2/060125_N5_o75

#directory containing fasta files of cds to be analyzed
HOG_CDS_DIR=${WD}/HOG_CDS

#list of HOG IDs from directory of fasta files
HOG_LIST=${HOG_CDS_DIR}/redo_hogs_max_rsc.txt

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
	java -jar -Xmx47G /home/crunnel2/bin/macse_v2.07.jar \
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