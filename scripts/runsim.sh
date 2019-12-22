#!/bin/bash

set -e

THREADS=1 # TODO user param

# Example
# ./scripts/runsim.sh -r /storage/resources/dbase/human/hg19/hg19.fa -p /storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB.bed -m /storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB-1.9.json -o del/del -x "--scale-outliers -c 7 --numcopies 100 --numreads 100000 --readlen 36 --region chr22:1-51304566"

SIMPATH=/storage/mgymrek/workspace/tulip/src/tulip

# Usage 
show_help () {
  cat <<ENDUSAGE
Usage:
  $(basename "$0") \\
     -r ref.fa \\
     -p peaks.bed \\
     -m model.json \\
     -o outprefix \\
     [-x extra params]

Script to run chips and all typical downstream analyses

Input: simulation parameters, reference genome, peaks, and output prefix.
Runs: (1) simulation, (2) alignment, (3) TDF conversion, (4) peak calling, (5) Output metrics

Also computes run time and memory usage of simulation.

ENDUSAGE
  exit 2
}

# Log
log() {
    BASE=$(basename "$0")
    echo "[$BASE]: $1" >&2
}

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# Parse options
OPTIND=1
while getopts "r:p:m:o:x:h?v" opt; do
    case "$opt" in
	h|\?)
            show_help
            exit 0
            ;;
	v) verbose=1
            ;;
	r) ref_fasta=$OPTARG
	    ;;
	p) peaks_file=$OPTARG
	    ;;
	m) model_file=$OPTARG
	    ;;
	o) output_prefix=$OPTARG
	    ;;
	x) extra_params=$OPTARG
	    ;;
    esac
done

shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift

# Check input options
if [[ -z "$ref_fasta" ]]; then
    die "No reference fasta specified (-r)"
fi

if [[ -z "$peaks_file" ]]; then
    die "No peaks file specified (-p)"
fi

if [[ -z "$model_file" ]]; then
    die "No model file specified (-m)"
fi

if [[ -z "$output_prefix" ]]; then
    die "No output prefix specified (-o)"
fi

log "Running with params: "
log "  verbose=$verbose"
log "  ref_fasta=${ref_fasta}"
log "  peaks_file=${peaks_file}"
log "  model_file=${model_file}"
log "  output_prefix=${output_prefix}"
log "  extra_params=${extra_params}"
log "  Left over unused options: $@"
log 

# Make sure output directory exists
mkdir -p $(dirname $output_prefix)

############### Step 1: Simulate ################
# Run simulator and capture timing output
# TODO get also memusg
log "  Running $SIMPATH"
{ time $SIMPATH simreads -p ${peaks_file} -t bed --model ${model_file} -f ${ref_fasta} -o ${output_prefix} ${extra_params} ; } 2> ${output_prefix}_time.log
time_real=$(cat ${output_prefix}_time.log | grep -w real | awk '{print $2}')
time_user=$(cat ${output_prefix}_time.log | grep -w user | awk '{print $2}')
time_sys=$(cat ${output_prefix}_time.log | grep -w sys | awk '{print $2}')

log "  Running $SIMPATH - WCE"
{ time $SIMPATH simreads -t wce --model ${model_file} -f ${ref_fasta} -o ${output_prefix}_wce ${extra_params} ; } 2> ${output_prefix}_timewce.log

############### Step 2: Alignment ################
# Run BWA-MEM. Different command based on paired vs. single
log "  Running bwa mem"
if [[ ${extra_params} == *"--paired"* ]]
then
    sh -c "bwa mem -t ${THREADS} ${ref_fasta} ${output_prefix}_1.fastq ${output_prefix}_2.fastq | samtools view -bS - > ${output_prefix}.bam " 
    sh -c "bwa mem -t ${THREADS} ${ref_fasta} ${output_prefix}_wce_1.fastq ${output_prefix}_wce_2.fastq | samtools view -bS - > ${output_prefix}_wce.bam " 
else
    sh -c "bwa mem -t ${THREADS} ${ref_fasta} ${output_prefix}.fastq | samtools view -bS - > ${output_prefix}.bam " 
    sh -c "bwa mem -t ${THREADS} ${ref_fasta} ${output_prefix}_wce.fastq | samtools view -bS - > ${output_prefix}_wce.bam " 
fi

# Sort and index the file
log "  Sorting and indexing the BAM"
samtools sort --threads $THREADS -T ${output_prefix} ${output_prefix}.bam > ${output_prefix}.sorted.bam 
samtools index ${output_prefix}.sorted.bam 

samtools sort --threads $THREADS -T ${output_prefix}_wce ${output_prefix}_wce.bam > ${output_prefix}_wce.sorted.bam 
samtools index ${output_prefix}_wce.sorted.bam 

# Mark duplicates
log "  Marking duplicates"
java -jar -Xmx12G -Djava.io.tmpdir=$(dirname ${output_prefix}) $PICARD MarkDuplicates \
    VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING \
    I=${output_prefix}.sorted.bam M=${output_prefix}.metrics O=${output_prefix}.flagged.bam 
java -jar -Xmx12G -Djava.io.tmpdir=$(dirname ${output_prefix}) $PICARD MarkDuplicates \
    VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING \
    I=${output_prefix}_wce.sorted.bam M=${output_prefix}_wce.metrics O=${output_prefix}_wce.flagged.bam 

############### Step 3: Make a TDF ################
log "  Make TDF"
igvtools count -z 5 -w 25 -e 0 ${output_prefix}.flagged.bam ${output_prefix}.tdf ${ref_fasta} 
igvtools count -z 5 -w 25 -e 0 ${output_prefix}_wce.flagged.bam ${output_prefix}_wce.tdf ${ref_fasta} 

############### Step 4: Peak calling ################
log "  Calling peaks"
macs2 callpeak --tempdir $(dirname ${output_prefix}) -t ${output_prefix}.flagged.bam -c ${output_prefix}_wce.flagged.bam --name ${output_prefix}_narrow  --nomodel
macs2 callpeak --broad --tempdir $(dirname ${output_prefix}) -t ${output_prefix}.flagged.bam -c ${output_prefix}_wce.flagged.bam --name ${output_prefix}_broad  --nomodel

############### Step 5: Compute metrics ################
numpeaks_orig=$(wc -l ${peaks_file})
numpeaks_narrow=$(wc -l ${output_prefix}_narrow_peaks.narrowPeak)
numpeaks_broad=$(wc -l ${output_prefix}_broad_peaks.broadPeak)

# SPOT-ORIG: perc of reads in original peaks - TODO
# SPOT-NARROW: perc of reads in narrow peaks - TODO
# SPOT-BROAD: perc of reacs in board peaks - TODO
# FPR-NARROW - TODO
# FPR-BROAD - TODO

# RECALL-NARROW: recall using narrow peaks
overlap_narrow=$(intersectBed -a ${peaks_file} -b ${output_prefix}_narrow_peaks.narrowPeak -u | wc -l)
recall_narrow=$(echo ${overlap_narrow}/${numpeaks_orig} | bc -l)
# RECALL-BROAD: recall using broad peaks
overlap_broad=$(intersectBed -a ${peaks_file} -b ${output_prefix}_broad_peaks.broadPeak -u | wc -l)
recall_broad=$(echo ${overlap_broad}/${numpeaks_orig} | bc -l)
# PEAKSIZE-NARROW: median peak size narrow
peaksize_narrow=$(cat ${output_prefix}_narrow_peaks.narrowPeak | awk '{print $3-$2}' | datamash median 1)
# PEAKSIZE-BROAD: median peak size broad
peaksize_broad=$(cat ${output_prefix}_broad_peaks.broadPeak | awk '{print $3-$2}' | datamash median 1)

###### Output results #####
echo "time-real,${time_real}" > ${output_prefix}_results.tab
echo "time-user,${time_user}" >> ${output_prefix}_results.tab
echo "time-sys,${time_sys}" >> ${output_prefix}_results.tab
echo "recall-narrow,${recall_narrow}" >> ${output_prefix}_results.tab
echo "recall-broad,${recall_broad}" >> ${output_prefix}_results.tab
echo "peaksize-narrow,${peaksize_narrow}" >> ${output_prefix}_results.tab
echo "peaksize-broad,${peaksize_broad}" >> ${output_prefix}_results.tab

