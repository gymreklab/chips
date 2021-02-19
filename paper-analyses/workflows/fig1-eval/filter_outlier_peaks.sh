#!/bin/bash

PEAKFILE=$1

# Set all peaks with scores >2x median score to have binding prob 1. Recommended with real data.
# TODO based on --scale-outliers

median=$(cat $PEAKFILE | datamash median 7)
#cat $PEAKFILE | awk -v"median=$median" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" (($7>2*median)?(100):($7*100/(2*median)))}'
cat $PEAKFILE

#q3=$(cat $PEAKFILE | datamash q3 7)
#iqr=$(cat $PEAKFILE | datamash iqr 7)
#mean=$(cat $PEAKFILE | datamash mean 7)
#cat $PEAKFILE | awk -v"q3=$q3" -v"iqr=$iqr" '($7<(q3+1.5*iqr))'
#max=$(cat $PEAKFILE | datamash max 7)
#cat $PEAKFILE | awk -v"max=$max" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7/(max)*100}'
