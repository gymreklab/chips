#!/bin/bash

OUTDIR=$1
SCALES=$2

REFFA=${OUTDIR}/chipmunk_spikein_1_1.fa

cat ${OUTDIR}/hg19.chr19_random_peaks.bed | awk '{print "TARGET-"$0}' | cut -f 1-3 > ${OUTDIR}/hg19.chr19_random_peaks_mod.bed

f=1
TAGDIRS=""
for scale in ${SCALES} 1
do
    mkdir -p ${OUTDIR}/chipmunk_spikein_${scale}_${f}
    makeTagDirectory ${OUTDIR}/chipmunk_spikein_${scale}_${f} ${OUTDIR}/chipmunk_spikein_${scale}_${f}.sorted.bam
    TAGDIRS="${TAGDIRS} ${OUTDIR}/chipmunk_spikein_${scale}_${f} "
done

annotatePeaks.pl ${OUTDIR}/hg19.chr19_random_peaks_mod.bed ${REFFA} \
    -size 1000 -hist 1 -d ${TAGDIRS} > ${OUTDIR}/spikein_histograms.txt
