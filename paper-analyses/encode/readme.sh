#!/bin/bash

OUTDIR=/storage/mgymrek/chipmunk/encode

# Run process encode examples
#thresh=5 # for HM
thresh=100 # for TF
while IFS='' read -r line || [[ -n "$line" ]]; do
    bamurl=$(echo $line | cut -f 4 -d',')
    bedurl=$(echo $line | cut -f 5 -d',')
    ct=$(echo $line | cut -f 1 -d',')
    f=$(echo $line | cut -f 2 -d',')
    bamacc=$(echo $bamurl | cut -d'/' -f 5)
    bedacc=$(echo $bedurl | cut -d'/' -f 5)
    factor=${ct}_${f}_${bamacc}_${bedacc}
#    echo ./process_encode.sh ${bamurl} ${bedurl} ${OUTDIR} ${factor} Single ${thresh}
    echo ./process_encode.sh ${bamurl} ${bedurl} ${OUTDIR} ${factor} Paired ${thresh}
done < encode_paired_example_datasets.csv | xargs -n1 -I% -P4 sh -c "%"
#encode_paired_example_datasets.csv #| xargs -n1 -I% -P4 sh -c "%"
#encode_datasets_K562_GM12878_clean_HM.csv
#encode_H3K27ac_reps.csv | xargs -n1 -I% -P4 sh -c "%"  
#< encode_datasets_K562_GM12878_clean_HM.csv | grep process | xargs -n1 -I% -P4 sh -c "%"
#< encode_paired_example_datasets.csv #| xargs -n1 -I% -P4 sh -c "%"

