#!/bin/bash

# Example:
# ./test_simlearn_consistency.sh 10 1 0.2
# e.g. Test that if we simulate s=0.2, we get out s=0.2

nc=$1
score=$2
spot=$3

# Generate random peaks
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from hg19.chromInfo" | grep -w chr19 > hg19.chr19.genome
bedtools random -l 200 -n 1000 -seed 100 -g hg19.chr19.genome | \
    awk -v"score=$score" '{print $1 "\t" $2 "\t" $3 "\t" score}' | \
    bedtools sort -i stdin > example_peaks.bed

# Simulate reads
./src/chipmunk simreads \
    -f /storage/resources/dbase/human/hg19/chromFa/chr19.fa \
    -p example_peaks.bed -t bed -c 4 \
    --gamma-frag 10,20 --spot $spot --frac 0.05 --pcr_rate 1.0 \
    --numreads 100000 --numcopies ${nc} --noscale --recomputeF \
    -o test.${nc} --readlen 100
bwa mem -t 10 /storage/resources/dbase/human/hg19/chromFa/chr19.fa \
    test.${nc}.fastq | samtools view -bS - > test.${nc}.bam
samtools sort -o test.${nc}.sorted.bam test.${nc}.bam
samtools index test.${nc}.sorted.bam

# Learn 
./src/chipmunk learn \
    -p example_peaks.bed -t bed -c 4 --thres 0 \
    -b test.${nc}.sorted.bam --noscale \
    -o test.${nc}
