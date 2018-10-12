#!/bin/bash

BAM=/storage/mlamkin/projects/encode_data/GM12878/ENCSR014YCR/released/hg19/alignments/bam/rep2/F1796-ENCFF037RQQ.sorted.flagged.bam
REFFA=/storage/resources/dbase/human/hg19/hg19.fa
OUTPREFIX=/storage/mgymrek/asimon-test
PEAKS=/storage/mlamkin/projects/encode_data/GM12878/ENCSR014YCR/released/hg19/alignments/bam/rep2/tags/peaks_lessthan404.txt
TYPE=homer
COLUMN=5

# spot: fraction of reads in peaks
# frac: fraction of genome that is bound

./src/asimon simreads \
    -p ${PEAKS} \
    -t ${TYPE} \
    -f ${REFFA} \
    -c ${COLUMN} \
    -o ${OUTPREFIX} \
    --numcopies 10000 \
    --numreads 100000 \
    --readlen 100 \
    --gamma-frag 5.663,44.08 \
    --spot 0.110237 --frac 0.00416413 \
    --region chr19:13628710-15628710 \
    --binsize 200000 \
    --paired

bowtie2 -x $(echo $REFFA | cut -d'.' -f 1) -1 ${OUTPREFIX}/reads_1.fastq -2 ${OUTPREFIX}/reads_2.fastq > ${OUTPREFIX}/reads.sam
samtools view -bS ${OUTPREFIX}/reads.sam > ${OUTPREFIX}/reads.bam
samtools sort ${OUTPREFIX}/reads.bam -o ${OUTPREFIX}/reads.sorted.bam
samtools index ${OUTPREFIX}/reads.sorted.bam
igvtools count -z 5 -w 25 -e 0 ${OUTPREFIX}/reads.sorted.bam  ${OUTPREFIX}/reads.sorted.tdf $REFFA

