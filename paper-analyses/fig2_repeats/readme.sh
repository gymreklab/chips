#!/bin/bash

source params.sh

# Make reference
./make_ref.sh

# Run simulation for various read lengths, paired/single
BAMS=""
TAGDIRS=""
for rl in 36 51 100 150
do
    #nreadsS=$(echo "${NUMTOTAL}/${rl}" | bc -l | cut -f 1 -d'.')
    #nreadsP=$(echo "${nreadsS}/2" | bc -l | cut -f 1 -d'.')
    nreadsS=100000
    nreadsP=100000
    ./run_repsim.sh ${rl} ${nreadsP} ${rl}.paired --paired
    ./run_repsim.sh ${rl} ${nreadsS} ${rl}.single
    BAMS="${BAMS} ${OUTDIR}/${rl}.single.sorted.bam ${OUTDIR}/${rl}.paired.sorted.bam "
    TAGDIRS="${TAGDIRS} ${OUTDIR}/${rl}.single/ ${OUTDIR}/${rl}.paired/"
done #| xargs -I% -P3 -n1 sh -c "%"

# Get coverage per bin + repeat lengths
multiBamCov -bams ${BAMS} -bed ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed | \
    intersectBed -b stdin -a ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed -wa -wb | \
    cut -f 6-10 --complement > ${OUTDIR}/repeat_cov_byrl.bed

# Get composite plots - all repeats
annotatePeaks.pl ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed ${CHROMFA} \
    -size 1000 -hist 1 \
    -d ${TAGDIRS} > ${OUTDIR}/repeat_composite.txt

# Get composite plots - repeats >40bp total
for thresh in 20 40 60 80 100
do
    cat ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed | awk -v"thresh=$thresh" '($5>=thresh)' > ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.gt${thresh}.bed
    annotatePeaks.pl ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.gt${thresh}.bed ${REFFA} \
    -size 1000 -hist 1 \
    -d ${TAGDIRS} > ${OUTDIR}/repeat_composite.gt${thresh}.txt
done
