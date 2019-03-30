#!/bin/bash

infile=$1

echo "| CellType|Factor|ModelFile |Model-k|Model-theta|Model-frac|Model-spot|Model-pcr| ENCODE_BAM|ENCODE_BED|"
echo "|---------|------|----------|----------|-------|----------|---------|---------|-----------|----------|"

cat $infile | grep -v CellType | \
    awk -F',' '{print $1 "|" $2 "|[json]("$NF ")|" $5"|"$6"|"$7"|"$8"|"$9"|"$3"|"$4}' | \
    sed 's/s3:\/\//http:\/\/s3.amazonaws.com\//' 
