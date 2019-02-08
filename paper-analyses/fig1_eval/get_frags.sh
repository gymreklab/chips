#!/bin/bash

source params.sh

factor=$1
threshold=200

BAM=$(ls /storage/mlamkin/projects/encode_data_learn/${factor}/*.bam)
BED=$(ls /storage/mlamkin/projects/encode_data_learn/${factor}/*.bed)

../../experimental/chip_sim_FINAL.py \
    $BAM \
    $BED \
    $threshold \
    ${OUTDIR}
