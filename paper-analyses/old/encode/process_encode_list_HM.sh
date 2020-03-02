#!/bin/bash

while IFS='' read -r line || [[ -n "$line" ]]; do
    # Check that not eGFP
    f=$(echo $line | cut -f 2 -d',')
    if [[ $f == *"eGFP-"* ]]; then
	continue
    fi
    # Check HM
    if [[ $f != *"H3"* ]]; then
	continue
    fi
    # Check cell type
    ct=$(echo $line | cut -f 1 -d',')
    if [[ $ct != "GM12878" ]] && [[ $ct != "K562" ]]; then
	continue
    fi
    # Check if bwa, no mark dup in header
    bamurl=$(echo $line | cut -f 4 -d',')
    header=$(samtools view -H ${bamurl})
    if [[ $header == *"MarkDuplicates"* ]]; then
	continue
    fi
    if [[ $header != *"bwa"* ]]; then
	continue
    fi
    echo $line
done < /storage/mlamkin/projects/encode_data_learn/datasets/learn_files/analyzed_encode_data.csv > encode_datasets_K562_GM12878_clean_HM.csv
