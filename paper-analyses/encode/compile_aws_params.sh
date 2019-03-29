#!/bin/bash

MDIR=/storage/mgymrek/chipmunk/encode/models-1.9

# Sync json results
aws s3 ls s3://chipmunk-encode-models/ | grep json | awk '{print $4}' | \
    grep "1\.9" | \
    xargs -n1 -P1 -I% sh -c "aws s3 cp s3://chipmunk-encode-models/% ${MDIR}/"

# Get local snorlax results 
find /storage/mgymrek/chipmunk/encode/ | grep json | grep "1\.9" | \
    xargs -n1 -P1 -I% sh -c "cp % ${MDIR}/"

echo "CellType,Factor,ENCODE_BAM,ENCODE_BED,Model-k,Model-theta,Model-frac,Model-spot,Model-pcr,ModelFile"
for f in $(ls $MDIR/*.json)
do
    ./print_model.py $f
done | sort -k1,1 -k2,2
