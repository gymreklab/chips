#!/bin/bash

set -e

# Run tests for varying read numbers for TF/HM

TFPEAKS=/storage/pandaman/project/ChIPs-experiments/data/TF-SP1-sim.bed
HMPEAKS=/storage/pandaman/project/ChIPs-experiments/data/HM-H3K27ac-sim.bed
TFMODEL=simplified_GM12878_CTCF.json
HMMODEL=simplified_GM12878_H3K27ac.json

for layout in single paired
do
    for style in "TF" "HM"
    do
	if [ x"$style" == x"TF" ]; then
	    MODELFILE=$TFMODEL
	    PEAKFILE=$TFPEAKS
	fi
	if [ x"$style" == x"HM" ]; then
	    MODELFILE=$HMMODEL
	    PEAKFILE=$HMPEAKS
	fi
	snakemake \
	    --config MODELFILE=$MODELFILE \
	    PEAKFILE=$PEAKFILE \
	    OUTPREFIX=${style}-${layout} \
	    LAYOUT=${layout} #--delete-all-output
    done
done
    
