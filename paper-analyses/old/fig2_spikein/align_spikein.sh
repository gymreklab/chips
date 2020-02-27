#!/bin/bash

OUTDIR=$1
SCALES=$2

REFFA=${OUTDIR}/chipmunk_spikein_1_1.fa

f=1
for scale in ${SCALES} 1
do
    bwa mem -t 10 ${REFFA} \
	${OUTDIR}/chipmunk_spikein_${scale}_${f}.fastq | \
	samtools view -bS - > ${OUTDIR}/chipmunk_spikein_${scale}_${f}.bam
    samtools sort -o ${OUTDIR}/chipmunk_spikein_${scale}_${f}.sorted.bam ${OUTDIR}/chipmunk_spikein_${scale}_${f}.bam
    samtools index ${OUTDIR}/chipmunk_spikein_${scale}_${f}.sorted.bam
done

