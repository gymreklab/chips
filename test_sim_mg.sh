#!/bin/bash

BAM=/storage/pandaman/project/asimon-data/ENCFF431YXJ.marked.bam
REFFA=/storage/resources/dbase/human/hg19/hg19.fa
OUTPREFIX=/storage/pandaman/project/asimon-test/reads
OUTDIR=/storage/pandaman/project/asimon-test/
PEAKS=/storage/pandaman/project/asimon-data/ENCFF878DLA.bed
TYPE=bed
SEQUENCER=HiSeq

# spot: fraction of reads in peaks
# frac: fraction of genome that is bound

#for nc in 1000 10000 100000 1000000
for nc in 10000
do
    ./src/chipmunk simreads \
	-p ${PEAKS} \
	-t ${TYPE} \
	-f ${REFFA} \
	-o ${OUTPREFIX} \
    -b ${BAM} \
	--numcopies ${nc} \
	--numreads 100000 \
	--readlen 100 \
	--gamma-frag 14.2142,15.9923 \
	--spot 0.00172121 --frac 8.43518e-05 \
	--region chr19:1-40000000 \
	--binsize 200000 \
    --thread 10 \
    --sequencer ${SEQUENCER}\
    --pcr_rate 0.8 \
	--paired
    
    bowtie2 -x $(echo $REFFA | cut -d'.' -f 1) -1 ${OUTDIR}/reads_1.fastq -2 ${OUTDIR}/reads_2.fastq > ${OUTDIR}/reads.${nc}.sam
    samtools view -bS ${OUTDIR}/reads.${nc}.sam > ${OUTDIR}/reads.${nc}.bam
    samtools sort ${OUTDIR}/reads.${nc}.bam -o ${OUTDIR}/reads.${nc}.sorted.bam
    samtools index ${OUTDIR}/reads.${nc}.sorted.bam
    igvtools count -z 5 -w 25 -e 0 ${OUTDIR}/reads.${nc}.sorted.bam  ${OUTDIR}/asimon.${nc}.sorted.tdf $REFFA
    rm ${OUTDIR}/reads*.fastq
done
