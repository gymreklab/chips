#!/bin/bash

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

# Run simulator and capture timing output
{ time $SIMPATH simreads -p ${peaks_file} -t bed --model ${model_file} -f ${ref_fasta} -o ${output_prefix} ${extra_params} ; } 2> ${output_prefix}_time.log
time_real=$(cat ${output_prefix}_time.log | grep -w real | awk '{print $2}')
time_user=$(cat ${output_prefix}_time.log | grep -w user | awk '{print $2}')
time_sys=$(cat ${output_prefix}_time.log | grep -w sys | awk '{print $2}')

###### Results #####
echo "Time: real=${time_real}; user=${time_user}; sys=${time_sys}"
    
