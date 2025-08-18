#!/bin/bash

for file in *
do
	export NAME=`basename $file ".cd-hit-est"`
	cd /media/will/mimic/transdecoder/clustered/TD-LO_complete_CR/$NAME/
	time TransDecoder.Predict -t /media/will/mimic/cd-hit/TD-P_complete_blast_pfam/$NAME.cd-hit-est \
		--retain_blastp_hits /media/will/mimic/transdecoder/clustered/blastp/$NAME.cd-hit-est.blastp.out
	mkdir /media/will/mimic/transdecoder/clustered/TD-P_complete_CR/$NAME/
	mv $NAME.* /media/will/mimic/transdecoder/clustered/TD-P_complete_CR/$NAME/
	mv pipeliner.* /media/will/mimic/transdecoder/clustered/TD-P_complete_CR/$NAME/
	# cd /media/will/mimic/cd-hit/TD-P_run_CR/
	# mv $NAME.cd-hit-est /media/will/mimic/cd-hit/TD-P_complete/
done
