#!/bin/bash

REFFA=/storage/resources/dbase/human/hg19/hg19.fa
PEAKS=testpeaks.tab

./src/asimon simreads \
    -p ${PEAKS} \
    -f ${REFFA} \
    -o test \
    --numcopies 10 \
    --numreads 1000 \
    --readlen 50 \
    --gamma-frag 200,0.5 \
    --spot 0.0035 --frac 0.0001 \
    --region chr20:45195673-45395941 \
    --binsize 10000
