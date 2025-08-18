#!/bin/bash

for file in *; do
export NAME=`basename $file .cd-hit-est`
echo $NAME
cd /media/will/mimic/transdecoder/clustered/blastp
time blastp -query /media/will/mimic/transdecoder/clustered/TD-LO_to_blastp/$NAME/$NAME.cd-hit-est.transdecoder_dir/longest_orfs.pep \
	-db /media/will/demogorgon/refdbs/uniprot/uniprot_sprot.fasta \
	-max_target_seqs 1 \
	-outfmt 6 \
	-evalue 1e-5 \
	-num_threads 36 \
	> $NAME.cd-hit-est.blastp.out
# need to cd back to directory where cd-hit files are here
mv $file /media/will/mimic/cd-hit/TD-P_run_CR/
done
