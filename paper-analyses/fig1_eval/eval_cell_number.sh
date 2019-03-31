#!/bin/bash

source params.sh

#binsize=5 # number of kb
#factor=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH

binsize=1
factor=GM12878_CTCF_ENCFF584BRF_ENCFF559IXF

# Setup windows on chr19
./make_windows.sh ${binsize}

#factor=GM12878_BACH1_ENCFF518TTP_ENCFF866OLZ

# Get genomic bins for sim/real
echo "Get encode bins"
#for f in $factor
#do
#    mkdir -p ${OUTDIR}/${f}/
#    BAM=${ENCDIR}/${f}/${f}.flagged.bam
#    bedtools multicov -bams $BAM -bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
#	${OUTDIR}/${f}/${f}.ENCODE.cov.${binsize}kb.bed
#    igvtools count ${BAM} ${OUTDIR}/${f}/${f}.encode.tdf ${REFFA}
#done

# Get overlap of bins with peak/no peak
BED=$(ls ${ENCDIR}/${factor}/*.bed | head -n 1) 
intersectBed \
    -a ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed \
    -b ${BED} -wa -u > ${OUTDIR}/${factor}/windows_${binsize}kb_peak.bed
intersectBed \
    -a ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed \
    -b ${BED} -wa -v > ${OUTDIR}/${factor}/windows_${binsize}kb_nopeak.bed

for nc in 1 5 10 25 50 100 1000
do
    echo $nc
    cat ${ENCDIR}/${factor}/*.json
#    ./run_factor.sh ${factor} ${nc} 51 269113 Single # info for h3k27ac
    ./run_factor.sh ${factor} ${nc} 36 285307 Single # info for CTCF. should be 28bp but doesn't work with bwa mem
#    ./run_factor.sh ${factor} ${nc} 101 665975 Paired # info for BACH1
    bedtools multicov -bams ${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam \
	-bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
	${OUTDIR}/${factor}/${factor}.${nc}.cov.${binsize}kb.bed
    ${CHIPMUNK} learn \
	-b ${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam \
	-p ${BED} -t bed -c 7 --thres 5 --scale-outliers \
	-o ${OUTDIR}/${factor}/${factor}.${nc}.json
done

