#!/bin/bash

BAM=/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/final-bams/F1796-HNFa3_S21.sorted.flagged.bam
PEAKS=/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/tagdir/F1796-HNFa3_S21/peaks.txt

./src/asimon learn \
    -b ${BAM} \
    -p ${PEAKS} \
    -o test
