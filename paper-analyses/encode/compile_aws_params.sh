#!/bin/bash

MDIR=/storage/mgymrek/chipmunk/encode/models

# Sync json results
cd $MDIR
aws s3 sync s3://chipmunk-encode-models .
cd -

for f in $(ls $MDIR/*.json)
do
    ./print_model.py $f
done
