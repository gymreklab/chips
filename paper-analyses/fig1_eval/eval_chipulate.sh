#!/bin/bash

source params.sh
binsize=5 # number of kb
factor=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH

# Setup windows on chr19
./make_windows.sh ${binsize}

BED=/storage/mgymrek/chipmunk/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH.bed
for nc in 1 2 5 10 100 1000 1000000
do
    echo $nc
    fq=/storage/yuq005/Chipmunk-eval/chipulate/run0327/reads_corrected/K27ac_chr19_${nc}.chip_reads.fastq.gz
    bwa mem -t 10 ${REFFA} ${fq} | samtools view -bS - > \
	${OUTDIR}/${factor}/${factor}.${nc}.chipulate.bam
    samtools sort -o ${OUTDIR}/${factor}/${factor}.${nc}.chipulate.sorted.bam ${OUTDIR}/${factor}/${factor}.${nc}.chipulate.bam
    samtools index ${OUTDIR}/${factor}/${factor}.${nc}.chipulate.sorted.bam
    java -jar $PICARD MarkDuplicates \
	I=${OUTDIR}/${factor}/${factor}.${nc}.chipulate.sorted.bam \
	O=${OUTDIR}/${factor}/${factor}.${nc}.chipulate.flagged.bam \
	M=${OUTDIR}/${factor}/${factor}.${nc}.chipulate.metrics \
	VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING
    samtools index ${OUTDIR}/${factor}/${factor}.${nc}.chipulate.flagged.bam
    igvtools count ${OUTDIR}/${factor}/${factor}.${nc}.chipulate.flagged.bam \
	${OUTDIR}/${factor}/${factor}.${nc}.chipulate.tdf ${REFFA}
    bedtools multicov -bams ${OUTDIR}/${factor}/${factor}.${nc}.chipulate.sorted.bam \
	-bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
	${OUTDIR}/${factor}/${factor}.${nc}.chipulate.cov.${binsize}kb.bed
done

# Get overlap of bins with peak/no peak
intersectBed \
    -a ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed \
    -b ${BED} -wa -u > ${OUTDIR}/${factor}/windows_${binsize}kb_peak.bed
intersectBed \
    -a ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed \
    -b ${BED} -wa -v > ${OUTDIR}/${factor}/windows_${binsize}kb_nopeak.bed
