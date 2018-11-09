#!/bin/bash

BAM=/storage/mlamkin/projects/encode_data/GM12878/ENCSR014YCR/released/hg19/alignments/bam/rep2/F1796-ENCFF037RQQ.sorted.flagged.bam
REFFA=/storage/resources/dbase/human/hg19/hg19.fa
OUTPREFIX=/storage/pandaman/project/asimon-test/reads
OUTDIR=/storage/pandaman/project/asimon-test/
PEAKS=/storage/mlamkin/projects/encode_data/GM12878/ENCSR014YCR/released/hg19/alignments/bam/rep2/tags/peaks_lessthan404.txt
TYPE=homer
COLUMN=5

# spot: fraction of reads in peaks
# frac: fraction of genome that is bound

for nc in 10 100 1000 1000 10000
do
    ./src/asimon simreads \
	-p ${PEAKS} \
	-t ${TYPE} \
	-f ${REFFA} \
	-c ${COLUMN} \
	-o ${OUTPREFIX} \
	--numcopies ${nc} \
	--numreads 100000 \
	--readlen 100 \
	--gamma-frag 5.663,44.08 \
	--spot 0.110237 --frac 0.00416413 \
	--region chr19:13628710-15628710 \
	--binsize 200000 \
	--paired
    
    bowtie2 -x $(echo $REFFA | cut -d'.' -f 1) -1 ${OUTDIR}/reads_1.fastq -2 ${OUTDIR}/reads_2.fastq > ${OUTDIR}/reads.${nc}.sam
    samtools view -bS ${OUTDIR}/reads.${nc}.sam > ${OUTDIR}/reads.${nc}.bam
    samtools sort ${OUTDIR}/reads.${nc}.bam -o ${OUTDIR}/reads.${nc}.sorted.bam
    samtools index ${OUTDIR}/reads.${nc}.sorted.bam
    igvtools count -z 5 -w 25 -e 0 ${OUTDIR}/reads.${nc}.sorted.bam  ${OUTDIR}/asimon.${nc}.sorted.tdf $REFFA
    rm ${OUTDIR}/reads*.fastq
done
