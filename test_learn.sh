#!/bin/bash

BAM=/storage/pandaman/project/asimon-data/ENCFF431YXJ.marked.bam
PEAKS=/storage/pandaman/project/asimon-data/ENCFF878DLA.bed
#BAM=/storage/mlamkin/projects/asimon-chip-sim/power-analysis/data/MCF-7/H3K27ac/ENCFF548YRU.flagged.bam
#PEAKS=/storage/mlamkin/projects/asimon-chip-sim/power-analysis/data/MCF-7/H3K27ac/ENCFF187RUK.bed
#BAM=/storage/mlamkin/projects/asimon-chip-sim/power-analysis/data/MCF-7/H3K9me3/ENCFF689XXY.flagged.bam
#PEAKS=/storage/mlamkin/projects/asimon-chip-sim/power-analysis/data/MCF-7/H3K9me3/ENCFF443NFW.bed
#BAM=/storage/mlamkin/projects/asimon-chip-sim/power-analysis/data/MCF-7/FOXA1/ENCFF287PVI.flagged.bam
#PEAKS=/storage/mlamkin/projects/asimon-chip-sim/power-analysis/data/MCF-7/FOXA1/ENCFF017WRM.bed
TYPE=bed
OUTDIR=/home/pandaman/projects/asimon/asimon-chip-sim/test/
THRES=250
COLUMN=7

./src/chipmunk learn \
    -p ${PEAKS} \
    -b ${BAM} \
    -t ${TYPE} \
    -o ${OUTDIR}\
    -c ${COLUMN}\
    --thres ${THRES}\
    --paired
