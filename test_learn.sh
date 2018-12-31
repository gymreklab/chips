#!/bin/bash

#BAM=/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/final-bams/F1796-HNFa3_S21.sorted.flagged.bam
#PEAKS=/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/tagdir/F1796-HNFa3_S21/peaks.txt #testpeaks.tab
BAM=/storage/pandaman/project/asimon-data/ENCFF431YXJ.marked.bam
PEAKS=/storage/pandaman/project/asimon-data/ENCFF878DLA.bed
TYPE=bed
OUTDIR=/home/pandaman/projects/asimon/asimon-chip-sim/test/

./src/asimon learn \
    -p ${PEAKS} \
    -b ${BAM} \
    -t ${TYPE} \
    -o ${OUTDIR}
