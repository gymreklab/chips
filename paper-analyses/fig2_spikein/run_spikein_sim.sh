#!/bin/bash

OUTDIR=$1
SCALES=$2

CHIPMUNK=chipmunk-1.9

f=1
for scale in ${SCALES} 1
do
    ${CHIPMUNK} simreads \
	-p ${OUTDIR}/chipmunk_spikein_${scale}_${f}.bed \
	-t bed -c 4 \
	-f ${OUTDIR}/chipmunk_spikein_1_1.fa \
	-o ${OUTDIR}/chipmunk_spikein_${scale}_${f} \
	--numcopies 100 --noscale --recomputeF \
	--binsize 1000000 \
	--numreads 1000000 \
	--model spikein_model.json
done
