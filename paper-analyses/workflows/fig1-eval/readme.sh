#!/bin/bash

# Make windows of chr21 of size 1kb, 5kb
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  | grep chr21 | head -n 1 > chr21.txt
bedtools makewindows -g chr21.txt -w 1000 > chr21_windows_1kb.bed
bedtools makewindows -g chr21.txt -w 5000 > chr21_windows_5kb.bed

PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH.bed
MODELFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH-1.9.json
OUTPREFIX=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH
snakemake \
    --config PEAKFILE=$PEAKFILE \
    MODELFILE=$MODELFILE \
    OUTPREFIX=$OUTPREFIX \
    LAYOUT=single \
    REF=/storage/resources/dbase/human/hg19/hg19.fa \
    REGION=chr22:1-51304566 \
    ENCDIR="" \
    C=7
