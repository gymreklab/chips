#!/bin/bash

source params.sh
factor=$1

GT=/storage/mgymrek/chipmunk_round2/encode/${factor}/${factor}.bed
numtotal=$(cat $GT | wc -l)

for numreads in $(echo $READNUMS | sed 's/,/ /g')
do
    f=/storage/mgymrek/chipmunk_round2/fig2_power/froms3/${factor}.${numreads}.${numreads}.broad_peaks.broadPeak
    numpeaks=$(intersectBed -a $GT -b $f -u | wc -l)
    echo $numreads $numpeaks $numtotal
done | awk -v"factor=$factor" '{print factor " " $0}' > ${OUTDIR}/${factor}.aws.broad.power.tab
