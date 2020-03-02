#!/bin/bash

source params.sh
factor=$1

GT=/storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed
numtotal=$(cat $GT | wc -l)

for numreads in $(echo $READNUMS | sed 's/,/ /g')
do
#    f=${OUTDIR}/${factor}.${numreads}.bed
    f=${OUTDIR}/${factor}.${numreads}.${numreads}_peaks.narrowPeak
    numpeaks=$(intersectBed -a $GT -b $f -u | wc -l)
    echo $numreads $numpeaks $numtotal
done | awk -v"factor=$factor" '{print factor " " $0}' > ${OUTDIR}/${factor}.power.tab
