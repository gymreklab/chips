#!/bin/bash

EXTRAARGS=$1 # e.g. --delete-all-output to clean

# Run all param workflows

#for usedir in "read-number" "read-length" "pcr" "frag-mean-length"
for usedir in "pcr"
do
    echo "##############" $usedir "#############"
    ./run.sh $usedir $EXTRAARGS
done
