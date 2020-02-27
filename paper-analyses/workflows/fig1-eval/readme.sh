#!/bin/bash

# Make windows of chr22 of size 1kb, 5kb
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  | grep chr22 | head -n 1 > chr22.txt
bedtools makewindows -g chr22.txt -w 1000 > chr22_windows_1kb.bed
bedtools makewindows -g chr22.txt -w 5000 > chr22_windows_5kb.bed

PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH.bed
MODELFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH-1.9.json
BAMFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH.flagged.bam
OUTPREFIX=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH

# Run simulations
snakemake \
    --config PEAKFILE=$PEAKFILE \
    MODELFILE=$MODELFILE \
    OUTPREFIX=$OUTPREFIX \
    BAMFILE=$BAMFILE \
    LAYOUT=single \
    REF=/storage/resources/dbase/human/hg19/hg19.fa \
    REGION=chr22:1-51304566 \
    ENCDIR="" \
    C=7
