#!/bin/bash

OUTDIR=$1
SCALES=$2

f=1
for scale in ${SCALES} 1
do
    # Get total # reads
    numreads=$(samtools view ${OUTDIR}/chipmunk_spikein_${scale}_${f}.sorted.bam | wc -l)
    # Get # reads mapped to drosophila
    numdros=$(samtools view ${OUTDIR}/chipmunk_spikein_${scale}_${f}.sorted.bam | grep SPIKE | wc -l)
    # Get # reads mapped to human
    numtar=$(samtools view ${OUTDIR}/chipmunk_spikein_${scale}_${f}.sorted.bam | grep TARGET | wc -l)
    echo ${scale} ${numreads} ${numdros} ${numtar}
done > ${OUTDIR}/chipmunk_spikein_readcounts.txt

