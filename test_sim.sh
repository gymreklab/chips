#!/bin/bash

BAM=/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/final-bams/F1796-HNFa3_S21.sorted.flagged.bam
REFFA=/storage/resources/dbase/human/hg19/hg19.fa
PEAKS=testpeaks.tab
TYPE=test
COLUMN=3

./src/asimon simreads \
    -p ${PEAKS} \
    -t ${TYPE} \
    -f ${REFFA} \
    -c ${COLUMN} \
    -o test \
    --numcopies 10 \
    --numreads 1000 \
    --readlen 50 \
    --gamma-frag 200,0.5 \
    --spot 0.0035 --frac 0.0001 \
    --region chr20:45195673-45395941 \
    --binsize 10000
