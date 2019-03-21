#!/bin/bash

source params.sh

# Make reference
./make_ref.sh

# Run simulation for various read lengths, paired/single
BAMS=""
TAGDIRS=""
for rl in 36 51 100 150
do
    nreadsS=1000000
    nreadsP=1000000
    ./run_repsim.sh ${rl} ${nreadsP} ${rl}.4.paired 4 " --paired --noscale"
    ./run_repsim.sh ${rl} ${nreadsS} ${rl}.4.single 4 " --noscale" 
    ./run_repsim.sh ${rl} ${nreadsP} ${rl}.5.paired 5 " --paired"
    ./run_repsim.sh ${rl} ${nreadsS} ${rl}.5.single 5 
    BAMS="${BAMS} ${OUTDIR}/${rl}.4.single.sorted.bam ${OUTDIR}/${rl}.4.paired.sorted.bam ${OUTDIR}/${rl}.5.single.sorted.bam ${OUTDIR}/${rl}.5.paired.sorted.bam"
    TAGDIRS="${TAGDIRS} ${OUTDIR}/${rl}.4.single/ ${OUTDIR}/${rl}.4.paired/ ${OUTDIR}/${rl}.5.single/ ${OUTDIR}/${rl}.5.paired/"
done #| xargs -I% -P3 -n1 sh -c "%"

# Get coverage per bin + repeat lengths
multiBamCov -bams ${BAMS} -bed ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed | \
    intersectBed -b stdin -a ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed -wa -wb | \
    cut -f 6-10 --complement > ${OUTDIR}/repeat_cov_byrl.bed

# Get composite plots long repeats
for thresh in 0 20 40 60 80 100
do
    cat ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed | awk -v"thresh=$thresh" '($5>=thresh)' | cut -f 1-3 > ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.gt${thresh}.bed
    annotatePeaks.pl ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.gt${thresh}.bed ${REFFA} \
	-size 500 -hist 1 \
	-d ${TAGDIRS} > ${OUTDIR}/repeat_composite.gt${thresh}.txt
done
