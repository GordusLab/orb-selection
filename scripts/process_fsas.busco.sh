#!/bin/bash

# run from /mimic/transdecoder/clustered/[directory with all TD-P outputs]

for dir in *
do
	export NAME=`echo $dir`
	cd /media/will/mimic/busco_260923/
	mkdir $NAME
	cd $NAME
	busco -i /media/will/mimic/transdecoder/clustered/TD-P_busco_run/$NAME/$NAME.cd-hit-est.transdecoder.pep \
		-l arachnida \
		-o $NAME.cd-hit-est.transdecoder.busco \
		-m prot \
		-c 60 \
		--datasets_version odb10
	mv /media/will/mimic/transdecoder/clustered/TD-P_busco_run/$NAME /media/will/mimic/transdecoder/clustered/TD-P_busco_complete/$NAME
done
