#!/bin/bash

source params.sh

factor=$1

BAM=$(ls /storage/mlamkin/projects/encode_data_learn/${factor}/*.bam)

python ../../experimental/view_pcr_dist.py $BAM
