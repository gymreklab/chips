#!/bin/bash

OUTDIR=/storage/mgymrek/chipmunk/encode

# Run process encode on 12 example factors
while IFS='' read -r line || [[ -n "$line" ]]; do
    bamurl=$(echo $line | cut -f 4 -d',')
    bedurl=$(echo $line | cut -f 5 -d',')
    ct=$(echo $line | cut -f 1 -d',')
    f=$(echo $line | cut -f 2 -d',')
    bamacc=$(echo $bamurl | cut -d'/' -f 5)
    bedacc=$(echo $bedurl | cut -d'/' -f 5)
    factor=${ct}_${f}_${bamacc}_${bedacc}
    echo ./process_encode.sh ${bamurl} ${bedurl} ${OUTDIR} ${factor} Both
done < encode_paired_example_datasets.csv #| xargs -n1 -I% -P4 sh -c "%"

