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
    if [ "x$LAYOUT" = "xpaired" ]
    then
	continue # only looking at single for chipulate right now
    fi

    echo $PEAKFILE

    # Get read info
    NUMREADS=$(samtools view $BAMFILE $CHROM | wc -l)
    NUMREADS=$((NUMREADS*2))
    READLEN=$(samtools view $BAMFILE | head -n 1000 | awk '{print length($10)}' | datamash max 1)
    READLEN=$(samtools view $BAMFILE | head -n 1000 | awk '{print length($10)}' | datamash max 1)

    # First without pcr
    snakemake $1 \
	--config PEAKFILE=$PEAKFILE \
	CHROM=$CHROM \
	MODELFILE=$OUTPREFIX.json \
	OUTPREFIX=$OUTPREFIX+"ischip" \
	BAMFILE=$BAMFILE \
	THRESH=$THRES \
	LAYOUT=$LAYOUT \
	REF=/storage/resources/dbase/human/hg19/hg19.fa \
	REFDIR=/storage/resources/dbase/human/hg19/chromFa \
	REGION=$REGION \
	ENCDIR="" \
	C=7 \
	NUMREADS=$NUMREADS \
	READLEN=$READLEN \
        METHOD=ischip PCR=0 HIGHBG=0
    # Now with PCR
    PCR=10
    snakemake $1 \
	--config PEAKFILE=$PEAKFILE \
	CHROM=$CHROM \
	MODELFILE=$OUTPREFIX.json \
	OUTPREFIX=$OUTPREFIX+"ischip-pcr${PCR}" \
	BAMFILE=$BAMFILE \
	THRESH=$THRES \
	LAYOUT=$LAYOUT \
	REF=/storage/resources/dbase/human/hg19/hg19.fa \
	REFDIR=/storage/resources/dbase/human/hg19/chromFa \
	REGION=$REGION \
	ENCDIR="" \
	C=7 \
	NUMREADS=$NUMREADS \
	READLEN=$READLEN \
	METHOD=ischip PCR=${PCR} HIGHBG=0
    
    # Now with high bg
    snakemake $1 \
	--config PEAKFILE=$PEAKFILE \
	CHROM=$CHROM \
	MODELFILE=$OUTPREFIX.json \
	OUTPREFIX=$OUTPREFIX+"ischip-HighBG" \
	BAMFILE=$BAMFILE \
	THRESH=$THRES \
	LAYOUT=$LAYOUT \
	REF=/storage/resources/dbase/human/hg19/hg19.fa \
        REFDIR=/storage/resources/dbase/human/hg19/chromFa \
	REGION=$REGION \
	ENCDIR="" \
	C=7 \
	NUMREADS=$NUMREADS \
	READLEN=$READLEN \
        METHOD=ischip PCR=0 HIGHBG=1

done < datasets.csv
