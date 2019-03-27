#!/bin/bash

source params.sh

factor=$1

./power-analysis.py \
    --bed /storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed \
    --model /storage/mgymrek/chipmunk/encode/${factor}/${factor}.json \
    --readnums ${READNUMS} \
    --threads 15 \
    --out ${OUTDIR}
