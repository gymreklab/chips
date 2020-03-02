#!/bin/bash

factor=$1
numreads=$2
optargs=$3

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# Set up
OUTDIR=/scratch/${factor}
mkdir -p ${OUTDIR} || die "Could not set up out dir"

# Get files needed
aws s3 cp s3://gymreklab/chipmunk-power/${factor}/${factor}.json ${OUTDIR}/${factor}.json || die "Could not get model file"
aws s3 cp s3://gymreklab/chipmunk-power/${factor}/${factor}.bed ${OUTDIR}/${factor}.bed || die "Could not get BED file"

# Get ref genome and index (if not there already)
if [ ! -f /scratch/hg19.fa ]; then
    aws s3 cp s3://gymreklab/resources/genomes/hg19/hg19.fa /scratch || die "Error copying ref genome"
    for suffix in amb ann bwt fai flat pac sa;
    do 
	aws s3 cp s3://gymreklab/resources/genomes/hg19/hg19.fa.${suffix} /scratch/ || die "Error copying bwa index"
    done
fi

# Run power analysis
power-analysis.py \
    --bed ${OUTDIR}/${factor}.bed \
    --model ${OUTDIR}/${factor}.json \
    --chipmunk-path chipmunk \
    --reffa /scratch/hg19.fa \
    --readnums ${numreads} \
    --out ${OUTDIR} --optargs "${optargs}"

# Upload results to s3.
aws s3 cp ${OUTDIR}/${factor}.${numreads}.${numreads}_peaks.narrowPeak s3://gymreklab/chipmunk-power/${factor}.${numreads}.${numreads}_peaks.narrowPeak
aws s3 cp ${OUTDIR}/${factor}.${numreads}.${numreads}.broad_peaks.broadPeak s3://gymreklab/chipmunk-power/${factor}.${numreads}.${numreads}.broad_peaks.broadPeak

# Tear down
rm -rf ${OUTDIR}
