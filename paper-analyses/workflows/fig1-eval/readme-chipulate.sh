#!/bin/bash

CHROM=chr22
REGION=chr22:1-51304566

while IFS= read -r line
do
    PEAKFILE=$(echo $line | cut -f 1 -d',')
    BAMFILE=$(echo $line | cut -f 2 -d',')
    OUTPREFIX=$(echo $line | cut -f 3 -d',')
#    THRES=$(echo $line | cut -f 4 -d',')
    THRES=0 # TODO check
    LAYOUT=$(echo $line | cut -f 5 -d',')
    if [ "x$PEAKFILE" = "xPEAKFILE" ]
    then
	continue
    fi
    echo $PEAKFILE

    # Run simulations
    NUMREADS=$(samtools view $BAMFILE $CHROM | wc -l)
    NUMREADS=$((NUMREADS*2))
    READLEN=$(samtools view $BAMFILE | head -n 1000 | awk '{print length($10)}' | datamash max 1)

    # Chipulate params
    NUMPEAKS=$(cat $PEAKFILE | grep -w ${CHROM} | wc -l)
    DEPTH=$(echo $NUMREADS/$NUMPEAKS | bc -l) # from chipulate.. equal to # reads/numpeaks

    READLEN=$(samtools view $BAMFILE | head -n 1000 | awk '{print length($10)}' | datamash max 1)

    if [ "x$LAYOUT" = "xpaired" ]
    then
	continue # only looking at single for chipulate right now
    fi    

    snakemake $1 \
	--config PEAKFILE=$PEAKFILE \
	CHROM=$CHROM \
	MODELFILE=$OUTPREFIX.json \
	OUTPREFIX=$OUTPREFIX+"chipulate" \
	BAMFILE=$BAMFILE \
	THRESH=$THRES \
	LAYOUT=$LAYOUT \
	REF=/storage/resources/dbase/human/hg19/hg19.fa \
	REGION=$REGION \
	ENCDIR="" \
	C=7 \
	NUMREADS=$NUMREADS \
	READLEN=$READLEN \
        METHOD=chipulate DEPTH=$DEPTH
    echo ""
done < datasets.csv
