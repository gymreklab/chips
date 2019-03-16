#!/bin/bash

BAMURL=$1
BEDURL=$2
FACTOR=$3
RTYPE=$4

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# Copy process encode script
aws s3 cp s3://gymreklab-awsbatch/process_encode.sh . || die "Couldn't get process_encode.sh"
chmod +x process_encode.sh || die "Couldn't chmod process_encode.sh"

# Clean out so we don't run out of space in case last job failed
mkdir -p /scratch/${FACTOR}/
rm -rf /scratch/${FACTOR}/*.

# Run process encode
echo "BAMURL: $BAMURL"
echo "BEDURL: $BEDURL"
echo "FACTOR: $FACTOR"
echo "RTYPE: $RTYPE"
./process_encode.sh ${BAMURL} ${BEDURL} /scratch ${FACTOR} ${RTYPE} || die "process_encode.sh failed"

# Upload results to s3
if [ "${RTYPE}" = "Single" ] || [ "${RTYPE}" = "Both" ]; then
    aws s3 cp /scratch/${FACTOR}/${FACTOR}.json s3://chipmunk-encode-models/${FACTOR}.json || die "Error writing results to s3"
fi

if [ "${RTYPE}" = "Paired" ] || [ "${RTYPE}" = "Both" ]; then
    aws s3 cp /scratch/${FACTOR}/${FACTOR}.paired.json s3://chipmunk-encode-models/${FACTOR}.paired.json || die "Error writing results to s3"
fi

# Remove the BAM file
rm -rf /scratch/${FACTOR}/*
