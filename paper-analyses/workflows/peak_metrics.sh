#!/bin/bash

# TODO: add recall as a function of peak score

TRUEPEAKS=$1
SIMPEAKS=$2
BAM=$3
REGION=$4 # Restrict analysis to this region

BLACKLIST=/storage/mgymrek/workspace/chips/paper-analyses/workflows/remove_regions_hg19.bed

numpeaks_sim=$(cat $SIMPEAKS | intersectBed -a stdin -b $BLACKLIST -v | wc -l)

# Recall. Only count appropriate region
if [[ -z "${REGION}" ]]
then
    numpeaks_true=$(cat $TRUEPEAKS | intersectBed -a stdin -b $BLACKLIST -v | wc -l)
else
    chrom=$(echo $REGION | cut -d':' -f 1)
    start=$(echo $REGION | cut -d':' -f 2 | cut -d'-' -f 1)
    end=$(echo $REGION | cut -d':' -f 2 | cut -d'-' -f 2)
    numpeaks_true=$(cat $TRUEPEAKS | intersectBed -a stdin -b $BLACKLIST -v | \
	awk -v"chrom=$chrom" -v "start=$start" -v "end=$end" '($1==chrom && $2>=start && $3<=end) {print $0}' | wc -l)
fi
overlap=$(intersectBed -a $TRUEPEAKS -b $SIMPEAKS -u | intersectBed -a stdin -b $BLACKLIST -v | wc -l)
recall=$(echo "${overlap}/${numpeaks_true}" | bc -l)
echo "recall,${recall}"

# False positive rate
wrong=$(intersectBed -a $SIMPEAKS -b $TRUEPEAKS -v | intersectBed -a stdin -b $BLACKLIST -v | wc -l)
fpr=$(echo "${wrong}/${numpeaks_sim}" | bc -l)
echo "fpr,${fpr}"

# Precision - TP/(TP+FP) = % of peaks that are true
sim_peaks_true=$(intersectBed -a $SIMPEAKS -b $TRUEPEAKS -u | intersectBed -a stdin -b $BLACKLIST -v | wc -l)
precision=$(echo ${sim_peaks_true}/${numpeaks_sim} | bc -l)
echo "precision,${precision}"

# Peak size
size=$(cat ${SIMPEAKS} | awk '{print $3-$2}' | datamash median 1)
echo "medpeaksize,${size}"

# SPOT (relative to true peaks)
numreads=$(samtools view -c $BAM)
inpeaks=$(intersectBed -abam $BAM -b ${TRUEPEAKS} | samtools view -c)
spot=$(echo "${inpeaks}/${numreads}" | bc -l)
echo "splot,${spot}"
