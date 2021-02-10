#!/bin/bash

############################################
# This snakemake workflow takes in an encode dataset (bed/bam) and:
# 1. Learns a model
# 2. Runs chips to simulate for various values of --numcopies
# 3. Aligns reads and flags dups
# 4. Counts reads in 1kb and 5kb windows for comparson with encode data

# Input datasets read from datasets.csv

# TODO: single vs. paired end mode for learn? add "layout" to datasets file
# TODO: set thresh automatically?
# TODO: add chipulate and isChIP for comparison
############################################

CHROM=chr22
REGION=chr22:1-51304566

while IFS= read -r line
do
    PEAKFILE=$(echo $line | cut -f 1 -d',')
    BAMFILE=$(echo $line | cut -f 2 -d',')
    OUTPREFIX=$(echo $line | cut -f 3 -d',')
    THRES=$(echo $line | cut -f 4 -d',')
    LAYOUT=$(echo $line | cut -f 5 -d',')
    if [ "x$PEAKFILE" = "xPEAKFILE" ]
    then
	continue
    fi
    echo $PEAKFILE

    # TODO  single or paired end mode? check ENCODE datasets above
    # Run simulations
    NUMREADS=$(samtools view $BAMFILE $CHROM | wc -l)
    NUMREADS=$((NUMREADS*2))
    READLEN=$(samtools view $BAMFILE | head -n 1000 | awk '{print length($10)}' | datamash max 1)
    snakemake \
	--config PEAKFILE=$PEAKFILE \
	CHROM=$CHROM \
	MODELFILE=$OUTPREFIX.json \
	OUTPREFIX=$OUTPREFIX \
	BAMFILE=$BAMFILE \
	THRESH=$THRES \
	LAYOUT=$LAYOUT \
	REF=/storage/resources/dbase/human/hg19/hg19.fa \
	REGION=$REGION \
	ENCDIR="" \
	C=7 \
	NUMREADS=$NUMREADS \
	READLEN=$READLEN \
        METHOD=chips
done < datasets.csv



