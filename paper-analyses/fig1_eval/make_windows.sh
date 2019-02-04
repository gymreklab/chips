#!/bin/bash

source params.sh

bedtools makewindows -g chr19.txt -w 1000 > ${OUTDIR}/windows/chr19_windows_1kb.bed
