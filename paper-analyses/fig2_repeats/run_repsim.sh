#!/bin/bash

READLEN=$1
NREADS=$2
PREFIX=$3
SCORECOL=$4
OPTARGS=$5

source params.sh

# Simulate reads
~/workspace/ChIPmunk/src/chipmunk simreads \
    -p ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed \
    -t bed -c ${SCORECOL} --recomputeF \
    -f ${CHROMFA} \
    -o ${OUTDIR}/${PREFIX} \
    --numcopies ${NC} \
    --gamma-frag ${GAMMA} \
    --spot ${SPOT} --pcr_rate 1 \
    --numreads ${NREADS} \
    --readlen ${READLEN} \
    --thread 10 \
    --region chr${CHROM}:1-249250621 $OPTARGS

# Map reads
if [[ $OPTARGS == *"--paired"* ]]; then
    bwa mem -t 10 ${REFFA} \
	${OUTDIR}/${PREFIX}_1.fastq \
	${OUTDIR}/${PREFIX}_2.fastq | \
	samtools view -bS - > ${OUTDIR}/${PREFIX}.bam
else
    bwa mem -t 10 ${REFFA} \
	${OUTDIR}/${PREFIX}.fastq | \
	samtools view -bS - > ${OUTDIR}/${PREFIX}.bam    
fi

# Sort and index
samtools sort -o ${OUTDIR}/${PREFIX}.sorted.bam ${OUTDIR}/${PREFIX}.bam
samtools index ${OUTDIR}/${PREFIX}.sorted.bam

# Convert to TDF
#igvtools count -w 1 ${OUTDIR}/${PREFIX}.sorted.bam ${OUTDIR}/${PREFIX}.tdf ${REFFA}

# Use Homer to make tag dir
mkdir -p ${OUTDIR}/${PREFIX}
makeTagDirectory ${OUTDIR}/${PREFIX} ${OUTDIR}/${PREFIX}.sorted.bam
