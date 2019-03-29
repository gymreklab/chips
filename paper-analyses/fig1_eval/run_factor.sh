#!/bin/bash
set -e
source params.sh

factor=$1
nc=$2
READLEN=$3
NREADS=$4
RTYPE=$5

BED=$(ls ${ENCDIR}/${factor}/*.bed | head -n 1)
BAM=$(ls ${ENCDIR}/${factor}/*.flagged.bam | head -n 1)
MODEL=$(ls ${ENCDIR}/${factor}/*.json | head -n 1)
if [ "$RTYPE" = "Paired" ]; then
    MODEL=$(ls ${ENCDIR}/${factor}/*.paired.json | head -n 1)
fi

mkdir -p ${OUTDIR}/${factor}

# Output params
echo "Readlen $READLEN"
echo "Nreads $NREADS"
echo "Rtype $RTYPE"

# Simulate reads
OPTARGS=""
if [ "$RTYPE" = "Paired" ]; then
    echo "Adding paired option"
    OPTARGS=" --paired"
fi
time $CHIPMUNK simreads \
    -p ${BED} \
    -t bed -c 7 --scale-outliers \
    -f ${REFFA} \
    -o ${OUTDIR}/${factor}/${factor}.${nc} \
    --model ${MODEL} \
    --numcopies ${nc} \
    --numreads ${NREADS} \
    --readlen ${READLEN} \
    --thread 10 \
    --region chr19:1-59128983 ${OPTARGS}

# Map reads
if [ "$RTYPE" = "Single" ]; then
    echo "Running single end"
    bwa mem -t 10 ${REFFA} \
	${OUTDIR}/${factor}/${factor}.${nc}.fastq | \
	samtools view -bS - > ${OUTDIR}/${factor}/${factor}.${nc}.bam
fi
if [ "$RTYPE" = "Paired" ]; then
    echo "Running paired end"
    bwa mem -t 10 ${REFFA} \
	${OUTDIR}/${factor}/${factor}.${nc}_1.fastq \
	${OUTDIR}/${factor}/${factor}.${nc}_2.fastq | \
	samtools view -bS - > ${OUTDIR}/${factor}/${factor}.${nc}.bam
fi

# Sort and index
samtools sort -o ${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam ${OUTDIR}/${factor}/${factor}.${nc}.bam
samtools index ${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam

# Mark duplicates
java -jar $PICARD MarkDuplicates \
    I=${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam \
    O=${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam \
    M=${OUTDIR}/${factor}/${factor}.${nc}.metrics \
    VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING
samtools index ${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam

# Convert to TDF
igvtools count ${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam \
    ${OUTDIR}/${factor}/${factor}.${nc}.tdf ${REFFA}
