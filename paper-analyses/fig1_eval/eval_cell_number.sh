#!/bin/bash

source params.sh
binsize=5 # number of kb

# Setup windows on chr19
./make_windows.sh ${binsize}

#factor=GM12878_BCLAF1_ENCFF800ZBG_ENCFF596IAY #GM12878_BACH1_ENCFF518TTP_ENCFF866OLZ
#factor2=GM12878_BCLAF1_ENCFF662GDH_ENCFF591DYQ
#factor=GM12878_MEF2A_ENCFF359HKJ_ENCFF980AUY
#factor2=GM12878_MEF2A_ENCFF012YAB_ENCFF951TLZ
factor=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH
factor2=GM12878_H3K27ac_ENCFF385RWJ_ENCFF849JNQ

# Get genomic bins for sim/real
echo "Get encode bins"
for f in $factor $factor2
do
    mkdir -p ${OUTDIR}/${f}/
    BAM=${ENCDIR}/${f}/${f}.flagged.bam
    bedtools multicov -bams $BAM -bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
	${OUTDIR}/${f}/${f}.ENCODE.cov.${binsize}kb.bed
    igvtools count ${BAM} ${OUTDIR}/${f}/${f}.encode.tdf ${REFFA}
done


for nc in 1 5 10 25 50 100 1000 10000
do
    echo $nc
    cat ${ENCDIR}/${factor}/*.json
    ./run_factor.sh ${factor} ${nc} 51 269113 Single #36 525661 Single # info for MEF2A #36 908621 Single # info for BCLAF1
    bedtools multicov -bams ${OUTDIR}/${factor}/${factor}.${nc}.sorted.bam \
	-bed ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed > \
	${OUTDIR}/${factor}/${factor}.${nc}.cov.${binsize}kb.bed
    BED=$(ls ${ENCDIR}/${factor}/*.bed | head -n 1) 
    ${CHIPMUNK} learn \
	-b ${OUTDIR}/${factor}/${factor}.${nc}.flagged.bam \
	-p ${BED} -t bed -c 7 --thres 5 \
	-o ${OUTDIR}/${factor}/${factor}.${nc}.json
done

