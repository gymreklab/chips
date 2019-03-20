#!/bin/bash

factor=GM12878_BCLAF1_ENCFF800ZBG_ENCFF596IAY
# Cat original model
cat /storage/mgymrek/chipmunk/encode/${factor}/${factor}.json
# Run learn on simulated dataset
chipmunk learn \
    -b /storage/mgymrek/chipmunk/fig1_eval/${factor}/${factor}.10000.flagged.bam \
    -p /storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed \
    -t bed -c 7 --thres 100 \
    -o test
