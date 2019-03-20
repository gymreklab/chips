#!/bin/bash

REFFA=$1
SPIKEFA=$2
OUTDIR=$3
SCALES=$4

# Randomly generate peaks from each ref genome
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from dm6.chromInfo" | grep chr4  > ${OUTDIR}/dm6.chr4.genome
bedtools random \
    -l 100 -n 20 -seed 100 -g ${OUTDIR}/dm6.chr4.genome | \
    awk '{print $1 "\t" $2 "\t" $3 "\t1"}' | \
    bedtools sort -i stdin > ${OUTDIR}/dm6.chr4_random_peaks.bed
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from hg19.chromInfo" | grep chr19 > ${OUTDIR}/hg19.chr19.genome
bedtools random -l 100 -n 10000 -seed 100 -g ${OUTDIR}/hg19.chr19.genome | \
    awk '{print $1 "\t" $2 "\t" $3 "\t1"}' | \
    bedtools sort -i stdin > ${OUTDIR}/hg19.chr19_random_peaks.bed

# Create spikein reference with no scaling
r=1; f=1
~/workspace/ChIPmunk/scripts/chipmunk-spike-in.sh \
    -t ${REFFA} \
    -s ${SPIKEFA} \
    -p ${OUTDIR}/hg19.chr19_random_peaks.bed \
    -q ${OUTDIR}/dm6.chr4_random_peaks.bed \
    -r ${r} -f ${f} \
    -o ${OUTDIR}/chipmunk_spikein_${r}_${f}

# Create spikein peaks for various scaling levels
f=1
for r in ${SCALES}
do
    ~/workspace/ChIPmunk/scripts/chipmunk-spike-in.sh \
	-t ${REFFA} \
	-s ${SPIKEFA} \
	-p ${OUTDIR}/hg19.chr19_random_peaks.bed \
	-q ${OUTDIR}/dm6.chr4_random_peaks.bed \
	-r ${r} -f ${f} -x \
	-o ${OUTDIR}/chipmunk_spikein_${r}_${f}
done
