#!/bin/bash

source params.sh

#binsize=1 # kb
#factor=GM12878_CTCF_ENCFF406XWF_ENCFF833FTF
#fqname="/storage/yuq005/Chipmunk-eval/chipulate/run0331_CTCF/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF_NUMCELLS.chip_reads.fastq"

binsize=5
factor=GM12878_H3K27ac_ENCFF385RWJ_ENCFF816AHV
fqname="/storage/yuq005/Chipmunk-eval/chipulate/run0331_K27ac/GM12878_H3K27ac_ENCFF385RWJ_ENCFF816AHV_NUMCELLS.chip_reads.fastq"

# Setup windows on chr19
./make_windows.sh ${binsize}

# Compare to chipmunk
nc=25
./run_factor.sh ${factor} ${nc} 36 521557 Single nobam
bedtools multicov -bams ${OUTDIR}/${factor}_nobam/${factor}_nobam.${nc}.flagged.bam \
    -bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
    ${OUTDIR}/${factor}_nobam/${factor}_nobam.${nc}.cov.${binsize}kb.bed

BED=/storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed
for nc in 1 2 5 10 100 1000 1000000
do
    echo $nc
    fq=$(echo $fqname | sed "s/NUMCELLS/${nc}/")
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
