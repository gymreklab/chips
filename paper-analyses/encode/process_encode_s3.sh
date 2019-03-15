#!/bin/bash

BAMURL=$1
BEDURL=$2
FACTOR=$3

# Copy process encode script
aws s3 cp s3://gymreklab-awsbatch/process_encode.sh .

# Run process encode
./process_encode.sh ${BAMURL} ${BEDURL} / ${FACTOR} Single

# Upload results to s3
aws s3 cp /${FACTOR}/${FACTOR}.json s3://chipmunk-encode-models/${FACTOR}.json
