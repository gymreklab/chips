#!/bin/bash

source params.sh

factor=$1
nc=$2

BED=$(ls ${ENCDIR}/${factor}/*.bed)
MODEL=$(ls ${ENCDIR}/${factor}/*.json | head -n 1)

mkdir ${OUTDIR}/${factor}
# Simulate reads
time $CHIPMUNK simreads \
    -p ${BED} \
    -t bed -c 7 \
    -f ${REFFA} \
    -o ${OUTDIR}/${factor}/${factor}.${nc} \
    --model ${MODEL} \
    --numcopies ${nc} \
    --numreads ${NREADS} \
    --readlen ${READLEN} \
    --region chr19:1-59128983 \
    --paired

# Map reads
bwa mem -t 10 ${REFFA} \
    ${OUTDIR}/${factor}/${factor}.${nc}_1.fastq \
    ${OUTDIR}/${factor}/${factor}.${nc}_2.fastq | \
    samtools view -bS - > ${OUTDIR}/${factor}.${factor}.${nc}.bam

# Sort and index
samtools sort -o ${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam ${OUTDIR}/${factor}/${factor}.${nc}.bam
samtools index ${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam

# Convert to TDF
igvtools count ${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam ${OUTDIR}/${factor}/${factor}.${nc}.tdf ${REFFA}
