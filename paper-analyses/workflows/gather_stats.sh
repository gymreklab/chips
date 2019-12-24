#!/bin/bash

# Input the three metrics files. Output key values in tab format
NARROW=$1
BROAD=$2
SIM=$3

echo "narrow-recall,narrow-fpr,narrow-medpeaksize,narrow-spot,broad-recall,broad-fpr,broad-medpeaksize,broad-spot,time_s,mem_MB" | \
    sed 's/,/\t/g'

output=$(cat $NARROW | grep recall | cut -d',' -f 2)
output=$output,$(cat $NARROW | grep fpr | cut -d',' -f 2)
output=$output,$(cat $NARROW | grep peaksize | cut -d',' -f 2)
output=$output,$(cat $NARROW | grep splot | cut -d',' -f 2)
output=$output,$(cat $BROAD | grep recall | cut -d',' -f 2)
output=$output,$(cat $BROAD | grep fpr | cut -d',' -f 2)
output=$output,$(cat $BROAD | grep peaksize | cut -d',' -f 2)
output=$output,$(cat $BROAD | grep splot | cut -d',' -f 2)
output=$output,$(cat $SIM | grep -v rss | awk -F"\t" '{print $1}')
output=$output,$(cat $SIM | grep -v rss | awk -F"\t" '{print $4}')
echo $output | sed 's/,/\t/g'
