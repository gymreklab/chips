#!/bin/bash

# Test command
# ./chipmunk-spike-in.sh -t /storage/resources/dbase/human/hg19/hg19.fa -s /storage/resources/dbase/drosophila/dm6/Sequence/WholeGenomeFasta/genome.fa -o /storage/mgymrek/chipmunk/fig2_spikein/spikein -f 1 -r 1 -x

# Initialize options
verbose=0
target_fasta="" # -t
spikein_fasta="" # -s
target_peaks="" # -p
spikein_peaks="" # -q
frac_spike=1 # -f
frac_target=1 # -r
output_prefix="" # -o
skip_ref=0

# Log
log() {
    BASE=$(basename "$0")
    echo "[$BASE]: $1" >&2
}

# Usage 
show_help () {
  cat <<ENDUSAGE
Usage:
  $(basename "$0") \\
     -t target.fa \\
     -s spikein.fa \\
     -p targetpeaks.bed \\
     -q spikepeaks.bed \\
     [-f spikefrac] \\
     [-x] \\
     -o outprefix

This script takes in the reference genome and peaks for the target genome and the spikein genome,
and outputs merged reference fasta and peak bed files to use for simulation. Assumes scores for 
peaks are in the 4th column of the peaks files and are between 0 and 1.

Required inputs:
-t <string>: Path to fasta reference file for the target genome
-s <string>: Path to fasta reference file for the spike in genome
-p <string>: Path to peaks for the target genome
-q <string>: Path to peaks for the spike in genome
-o <string>: Prefix to name output files

Optional arguments:
-f <float>: Scale peak scores for the spikein.
-r <float>: Scale peak scores for the target.
-x: Skip reference fasta creation, just merge the peak files

ENDUSAGE
  exit 2
}

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# Parse options
OPTIND=1
while getopts "t:s:p:q:f:r:o:xh?v" opt; do
    case "$opt" in
	h|\?)
            show_help
            exit 0
            ;;
	v) verbose=1
            ;;
	t) target_fasta=$OPTARG
	    ;;
	s) spikein_fasta=$OPTARG
	    ;;
	p) target_peaks=$OPTARG
	    ;;
	q) spikein_peaks=$OPTARG
	    ;;
	f) frac_spike=$OPTARG
	    ;;
	r) frac_target=$OPTARG
	    ;;
	o) output_prefix=$OPTARG
            ;;
	x) skip_ref=1
	    ;;
    esac
done

shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift

# Check input options
if [[ -z "$target_fasta" ]]; then
    die "No target fasta specified (-t)"
fi

if [[ -z "$spikein_fasta" ]]; then
    die "No spikein fasta specified (-s)"
fi

if [[ -z "$target_peaks" ]]; then
    die "No target peaks specified (-p)"
fi

if [[ -z "$spikein_peaks" ]]; then
    die "No spikein peaks specified (-q)"
fi

if [[ -z "$output_prefix" ]]; then
    die "No output prefix specified (-o)"
fi

log "Running with params: "
log "  verbose=$verbose"
log "  target_fasta=${target_fasta}"
log "  spikein_fasta=${spikein_fasta}"
log "  target_peaks=${target_peaks}"
log "  spikein_peaks=${spikein_peaks}"
log "  output_prefix=${output_prefix}"
log "  frac_spike=${frac_spike}"
log "  frac_target=${frac_target}"
log "  Left over unused options: $@"
log 

# Concatenate fastas
if [ "${skip_ref}" = 0 ]; then
    log "Making merged FASTA reference"
    cat ${target_fasta} | sed 's/^>/>TARGET-/' > ${output_prefix}.fa
    cat ${spikein_fasta} | sed 's/^>/>SPIKEIN-/' >> ${output_prefix}.fa
    samtools faidx ${output_prefix}.fa
    log "Output merged reference to ${output_prefix}.fa"
else
    log "Skipping FASTA merging"
fi

# Concatenate Peak files
log "Making merged peaks file"
cat ${target_peaks} | \
    awk -v"scale=${frac_target}" '{print "TARGET-"$1 "\t" $2 "\t" $3 "\t" $4*scale}' \
    > ${output_prefix}.bed
cat ${spikein_peaks} | \
    awk -v"scale=${frac_spike}" '{print "SPIKEIN-"$1 "\t" $2 "\t" $3 "\t" $4*scale}' \
    >> ${output_prefix}.bed
log "Output merged peaks to ${output_prefix}.bed"

log "Done preparing spikein files for input to chipmunk simreads!"
exit 0
