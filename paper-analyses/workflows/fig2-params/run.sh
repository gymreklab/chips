#!/bin/bash

set -e

USEDIR=$1 # run snakemake in this directory
EXTRAARGS=$2 # e.g. --delete-all-output 

# Run tests for varying read numbers for TF/HM

TFPEAKS=/storage/pandaman/project/ChIPs-experiments/data/TF-SP1-sim.bed
HMPEAKS=/storage/pandaman/project/ChIPs-experiments/data/HM-H3K27ac-sim.bed
TFMODEL=$(pwd)/simplified_TF.json
HMMODEL=$(pwd)/simplified_HM.json

cd $USEDIR

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
	    LAYOUT=${layout} ${EXTRAARGS}
    done
done
