#!/bin/bash


CHROM=chr22
REGION=chr22:1-51304566

while IFS= read -r line
do
    PEAKFILE=$(echo $line | cut -f 1 -d',')
    BAMFILE=$(echo $line | cut -f 2 -d',')
    OUTPREFIX=$(echo $line | cut -f 3 -d',')
    #THRES=$(echo $line | cut -f 4 -d',')
    THRES=0 # TODO check
    LAYOUT=$(echo $line | cut -f 5 -d',')
    if [ "x$PEAKFILE" = "xPEAKFILE" ]
    then
	continue
    fi
    echo $PEAKFILE

    # Run simulations
    NUMREADS=$(samtools view $BAMFILE $CHROM | wc -l)
    if [ "x$LAYOUT" = "xsingle" ]
    then
	NUMREADS=$((NUMREADS*2))
    fi

    READLEN=$(samtools view $BAMFILE | head -n 1000 | awk '{print length($10)}' | datamash max 1)
    snakemake $1 \
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
