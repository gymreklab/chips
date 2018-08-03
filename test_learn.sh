#!/bin/bash

BAM=/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/final-bams/F1796-HNFa3_S21.sorted.flagged.bam
PEAKS=testpeaks.tab

./src/asimon learn \
    -b ${BAM} \
    -p ${PEAKS} \
    -o test
