#!/bin/bash

source params.sh

REPBED=/storage/mgymrek/fromMartinOrth/GGAAmSats_hg19_tabdelimited.bed
cat ${REPBED} | \
    awk -v"ex=${PEXTEND}" '{print "chr"$1 "\t" int(($2+$3)/2)-ex "\t" int(($2+$3)/2)+ex "\t" 1 "\t" ($3-$2+1)}' | grep -w chr${CHROM} \
    > ${OUTDIR}/hg19.hipstr.chr${CHROM}.AAGG.bed
