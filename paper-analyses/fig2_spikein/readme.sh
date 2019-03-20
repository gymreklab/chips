#!/bin/bash

set -e

REFFA=/storage/resources/dbase/human/hg19/chromFa/chr19.fa
SPIKEFA=/storage/resources/dbase/drosophila/dm6/Sequence/Chromosomes/chr4.fa
OUTDIR=/storage/mgymrek/chipmunk/fig2_spikein/
SCALES="0.01 0.1 0.25 0.5 0.8"

# Get peaks with different scalings
echo "Generate peaks and reffa"
./generate_peaks.sh ${REFFA} ${SPIKEFA} ${OUTDIR} "${SCALES}"
bwa index ${OUTDIR}/chipmunk_spikein_1_1.fa

# Run simulation and alignment for each
echo "Run simulation"
./run_spikein_sim.sh ${OUTDIR} "${SCALES}"

# Align reads
echo "Align reads"
./align_spikein.sh ${OUTDIR} "${SCALES}"

# Get tag histograms
echo "Generate tag histograms"
./make_hists.sh ${OUTDIR} "${SCALES}"

# Get read counts
echo "Get read counts"
./get_read_counts.sh ${OUTDIR} "${SCALES}" 
