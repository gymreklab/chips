#!/bin/bash

source params.sh

factor=$1
nc=$2

BAM=$(ls /storage/mlamkin/projects/encode_data_learn/${factor}/*.bam)
BED=$(ls /storage/mlamkin/projects/encode_data_learn/${factor}/*.bed)

# Simulate reads
~/workspace/asimon-chip-sim/src/asimon simreads \
    -p ${BED} \
    -t bed -c 6 \
    -f ${REFFA} \
    -o ${OUTDIR}/${factor}.${nc} \
    --numcopies ${nc} \
    --gamma-frag 15.9201,14.9829 \
    --spot 0.201711 --frac 0.0425032 \
    --numreads ${NREADS} \
    --readlen ${READLEN} \
    --thread 1 \
    --region chr19:1-59128983 \
    --paired

# Map reads
bwa mem -t 10 ${REFFA} \
    ${OUTDIR}/${factor}.${nc}_1.fastq \
    ${OUTDIR}/${factor}.${nc}_2.fastq | \
    samtools view -bS - > ${OUTDIR}/${factor}.${nc}.bam

# Sort and index
samtools sort -o ${OUTDIR}/${factor}.${nc}.sorted.bam ${OUTDIR}/${factor}.${nc}.bam
samtools index ${OUTDIR}/${factor}.${nc}.sorted.bam

# Convert to TDF
igvtools count ${OUTDIR}/${factor}.${nc}.sorted.bam ${OUTDIR}/${factor}.${nc}.tdf ${REFFA}
igvtools count ${BAM} ${OUTDIR}/${factor}.ENCODE.tdf ${REFFA}

# Get genomic bins for sim/real
bedtools multicov -bams $BAM -bed ${OUTDIR}/windows/chr19_windows_1kb.bed > \
    ${OUTDIR}/${factor}.ENCODE.cov.1kb.bed
bedtools multicov -bams ${OUTDIR}/${factor}.${nc}.sorted.bam \
    -bed ${OUTDIR}/windows/chr19_windows_1kb.bed > \
    ${OUTDIR}/${factor}.${nc}.cov.1kb.bed
