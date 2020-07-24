#!/bin/bash

# Make windows of chr22 of size 1kb, 5kb
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  | grep chr22 | head -n 1 > chr22.txt
#bedtools makewindows -g chr22.txt -w 1000 > chr22_windows_1kb.bed
#bedtools makewindows -g chr22.txt -w 5000 > chr22_windows_5kb.bed

# Make windows of chr19 of size 1kb, 5kb
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  | grep chr19 | head -n 1 > chr19.txt
bedtools makewindows -g chr19.txt -w 1000 > chr19_windows_1kb.bed
bedtools makewindows -g chr19.txt -w 5000 > chr19_windows_5kb.bed


#PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH.bed
#BAMFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH/GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH.flagged.bam
#OUTPREFIX=GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH
#THRES=50
#READLEN=36


PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF.bed
BAMFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF.flagged.bam
MODELFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF/GM12878_CTCF_ENCFF406XWF_ENCFF833FTF-1.9.json
OUTPREFIX=GM12878_CTCF_ENCFF406XWF_ENCFF833FTF
THRES=5000
READLEN=36


#PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ/GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ.bed
#BAMFILE=/storage/mgymrek/chipmunk_round2/encode/GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ/GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ.flagged.bam
#OUTPREFIX=GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ
#THRES=20
#READLEN=36

#PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/K562_H3K4me3_ENCFF681JQI_ENCFF127XXD/K562_H3K4me3_ENCFF681JQI_ENCFF127XXD.bed
#BAMFILE=/storage/mgymrek/chipmunk_round2/encode/K562_H3K4me3_ENCFF681JQI_ENCFF127XXD/K562_H3K4me3_ENCFF681JQI_ENCFF127XXD.flagged.bam
#OUTPREFIX=K562_H3K4me3_ENCFF681JQI_ENCFF127XXD
#THRES=5
#READLEN=36

#PEAKFILE=/storage/mgymrek/chipmunk_round2/encode/K562_NFYA_ENCFF000YUR_ENCFF003WYE/K562_NFYA_ENCFF000YUR_ENCFF003WYE.bed
#BAMFILE=/storage/mgymrek/chipmunk_round2/encode/K562_NFYA_ENCFF000YUR_ENCFF003WYE/K562_NFYA_ENCFF000YUR_ENCFF003WYE.flagged.bam
#OUTPREFIX=K562_NFYA_ENCFF000YUR_ENCFF003WYE
#THRES=50
#READLEN=36

# Make model
#../../../src/chips learn \
#    -b $BAMFILE \
#    -p $PEAKFILE \
#    -t bed \
#    -c 7 \
#    -o $OUTPREFIX \
#    --scale-outliers \
#    --thres $THRES

# Run simulations
#NUMREADS=$(samtools view $BAMFILE chr22 | wc -l)
#NUMREADS=$((NUMREADS*2))
NUMREADS=521557
snakemake \
    --config PEAKFILE=$PEAKFILE \
    MODELFILE=$MODELFILE \
    OUTPREFIX=$OUTPREFIX \
    BAMFILE=$BAMFILE \
    LAYOUT=single \
    REF=/storage/resources/dbase/human/hg19/hg19.fa \
    REGION=chr19:1-59128983 \
    ENCDIR="" \
    C=7 \
    NUMREADS=$NUMREADS \
    READLEN=$READLEN
