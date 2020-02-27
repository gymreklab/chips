#!/bin/bash

OUTDIR=/storage/mgymrek/chipmunk/encode

# Run process encode examples
while IFS='' read -r line || [[ -n "$line" ]]; do
    ct=$(echo $line | cut -f 1 -d',')
    f=$(echo $line | cut -f 2 -d',')
    rtype=$(echo $line | cut -f 3 -d',')
    bamurl=$(echo $line | cut -f 4 -d',')
    bedurl=$(echo $line | cut -f 5 -d',')
    bamacc=$(echo $bamurl | cut -d'/' -f 5)
    bedacc=$(echo $bedurl | cut -d'/' -f 5)
    factor=${ct}_${f}_${bamacc}_${bedacc}
    if [[ "$factor" == *"H3"* ]]; then
	thresh=5 # For HMs
    else
	thresh=100 # For TFs
    fi
    if [[ "$rtype" == "Paired" ]]; then
	rtype=Both
    fi
    echo ./process_encode.sh ${bamurl} ${bedurl} ${OUTDIR} ${factor} ${rtype} ${thresh}
done < encode_examples_snorlax.csv #| xargs -n1 -I% -P4 sh -c "%" 
