#!/bin/bash

source params.sh

binsize=$1
factor=$2
readlen=$3
numreads=$4
rtype=$5

#binsize=5 # kb
#factor=GM12878_H3K27ac_ENCFF385RWJ_ENCFF816AHV
#readlen=51
#numreads=269113
#rtype=Single

#binsize=1 # kb
#factor=GM12878_CTCF_ENCFF406XWF_ENCFF833FTF
#readlen=36
#numreads=521557
#rtype=Single

mkdir -p ${OUTDIR}/${factor}/

# Setup windows on chr19
./make_windows.sh ${binsize}

# Get genomic bins for sim/real
echo "Get encode bins"
f=$factor
BAM=${ENCDIR}/${f}/${f}.flagged.bam
bedtools multicov -bams $BAM -bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
    ${OUTDIR}/${f}/${f}.ENCODE.cov.${binsize}kb.bed
igvtools count ${BAM} ${OUTDIR}/${f}/${f}.encode.tdf ${REFFA}

# Get overlap of bins with peak/no peak
BED=$(ls ${ENCDIR}/${factor}/*.bed | head -n 1) 
intersectBed \
    -a ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed \
    -b ${BED} -wa -u > ${OUTDIR}/${factor}/windows_${binsize}kb_peak.bed
intersectBed \
    -a ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed \
    -b ${BED} -wa -v > ${OUTDIR}/${factor}/windows_${binsize}kb_nopeak.bed

#for nc in 1 5 10 25 50 100 1000 5000
for nc in 100
do
    echo $nc
    cat ${ENCDIR}/${factor}/*.json
    ./run_factor.sh ${factor} ${nc} ${readlen} ${numreads} ${rtype}
    bedtools multicov -bams ${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam \
	-bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
	${OUTDIR}/${factor}/${factor}.${nc}.cov.${binsize}kb.bed
done

