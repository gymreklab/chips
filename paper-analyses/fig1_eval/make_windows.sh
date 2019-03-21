#!/bin/bash

source params.sh

binsize=$1

bedtools makewindows -g chr19.txt -w ${binsize}000 > ${OUTDIR}/windows/chr19_windows_${binsize}kb.bed
