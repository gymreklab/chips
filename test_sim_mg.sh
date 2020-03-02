#!/bin/bash

BAM=/storage/pandaman/project/asimon-data/ENCFF431YXJ.marked.bam
REFFA=/storage/resources/dbase/human/hg19/hg19.fa
OUTPREFIX=/storage/pandaman/project/asimon-test/reads
OUTDIR=/storage/pandaman/project/asimon-test/
PEAKS=/storage/pandaman/project/asimon-data/ENCFF878DLA.bed
#PEAKS=/storage/pandaman/project/asimon-data/test.bed
TYPE=bed
SEQUENCER=HiSeq

# spot: fraction of reads in peaks
# frac: fraction of genome that is bound

#for nc in 1000 10000 100000 1000000
for nc in 1000
do
    ./src/chips simreads \
	-p ${PEAKS} \
	-t ${TYPE} \
	-f ${REFFA} \
	-o ${OUTPREFIX} \
    -b ${BAM} \
	--numcopies ${nc} \
	--numreads 100000 \
	--readlen 100 \
	--gamma-frag 15.2896,15.0162 \
	--spot 0.245845 --frac 0.000179822 \
	--binsize 200000 \
    --thread 15 \
    --region chr19:1-40000000 \
    --sequencer ${SEQUENCER}\
    --pcr_rate 0.825219 \
    --seed 396587666 \
	--paired
    
    bowtie2 -x $(echo $REFFA | cut -d'.' -f 1) -1 ${OUTDIR}/reads_1.fastq -2 ${OUTDIR}/reads_2.fastq > ${OUTDIR}/reads.${nc}.sam
    samtools view -bS ${OUTDIR}/reads.${nc}.sam > ${OUTDIR}/reads.${nc}.bam
    samtools sort ${OUTDIR}/reads.${nc}.bam -o ${OUTDIR}/reads.${nc}.sorted.bam
    samtools index ${OUTDIR}/reads.${nc}.sorted.bam
    igvtools count -z 5 -w 25 -e 0 ${OUTDIR}/reads.${nc}.sorted.bam  ${OUTDIR}/asimon.${nc}.sorted.tdf $REFFA
    rm ${OUTDIR}/reads*.fastq
done
