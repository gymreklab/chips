#!/bin/bash

OUTDIR=/storage/mgymrek/chipmunk/encode

# Run process encode on factors
# TODO for loop to run
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo $line
    #./process_encode.sh ${bamurl} ${bedurl} ${OUTDIR} ${factor} Both
done < encode_paired_example_datasets.csv
