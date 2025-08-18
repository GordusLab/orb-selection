#!/bin/bash

# run from directory containing cd-hit files

for file in *; do
export NAME=`basename $file .cd-hit-est`
cd /media/will/mimic/transdecoder/clustered/TD-LO_complete_CR/
mkdir $NAME
cd /media/will/mimic/transdecoder/clustered/TD-LO_complete_CR/$NAME
TransDecoder.LongOrfs -t /media/will/mimic/cd-hit/TD-P_complete_blast_pfam/$NAME.cd-hit-est
# this script doesn't move the cd-hit files
done
