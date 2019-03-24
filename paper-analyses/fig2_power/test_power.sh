#!/bin/bash

source params.sh

factor=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH

READNUMS=1000000

./power-analysis.py \
    --bed /storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed \
    --model /storage/mgymrek/chipmunk/encode/${factor}/${factor}.json \
    --readnums ${READNUMS} \
    --out ${OUTDIR} --optargs " --region chr19:1-59128983" --debug
