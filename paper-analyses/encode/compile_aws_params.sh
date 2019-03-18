#!/bin/bash

MDIR=/storage/mgymrek/chipmunk/encode/models

# Sync json results
cd $MDIR
aws s3 sync s3://chipmunk-encode-models . > /dev/null 2>&1 
cd -

echo "CellType,Factor,ENCODE_BAM,ENCODE_BED,Model-k,Model-theta,Model-frac,Model-spot,Model-pcr,ModelFile"
for f in $(ls $MDIR/*.json)
do
    ./print_model.py $f
done | sort -k1,1 -k2,2
