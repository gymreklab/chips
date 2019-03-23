#!/bin/bash

# Test
factor=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH
./power-analysis.py \
    --bed /storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed \
    --model /storage/mgymrek/chipmunk/encode/${factor}/${factor}.json \
    --readnums 1000000,5000000,10000000,50000000 \
    --out /storage/mgymrek/chipmunk/fig2_power/
