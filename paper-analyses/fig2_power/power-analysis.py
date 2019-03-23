#!/usr/bin/env python3
"""
ChIPmunk power analysis

"""
# Author Michael Lamkin Last Edited 3/8/19
# Edited by Melissa Gymrek 3/22/19

import argparse
import sys, os
from subprocess import Popen, PIPE
import glob
from collections import defaultdict

CHIPMUNK="/home/mgymrek/workspace/ChIPmunk/src/chipmunk"
MODEL_K = 10
MODEL_THETA = 20
MODEL_PCR = 0.85
PARAM_NC = 100
PARAM_RL = 100
THREADS = 10
REFFA = '/storage/resources/dbase/human/hg19/hg19.fa'
RLIST="100000,500000,1000000,5000000,1e8"

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--bed", help="BED file with peaks", type=str, required=True)
    parser.add_argument("--model", help="JSON model file", type=str, required=True)
    parser.add_argument("--readnums", help="comma separated list of read nums", type=str, default=RLIST)
    parser.add_argument("--out", help="Output directory", type=str, required=True)

    args = parser.parse_args()
    OUTDIR = args.out
    read_list = [int(item) for item in args.readnums.split(",")]

    for numreads in read_list:
        bed_file = args.bed
        factor = os.path.basename(bed_file).strip(".bed")
        output_prefix = os.path.join(OUTDIR, factor+"."+str(numreads))
#                     '--region chr19:1-59128983 '
        chipmunk = ('%s simreads '
                    '-p %s '
                    '-t bed '
                    '-f %s '
                    '-c 7 '
                    '-o %s '
                    '--numcopies %d '
                    '--numreads %d '
                    '--readlen %d '
                    '--model %s '
                    '--gamma-frag %s,%s '
                    '--pcr_rate %f '
                    '--scale-outliers '
                    '--thread %d ')%(CHIPMUNK, bed_file, REFFA, output_prefix, PARAM_NC, numreads, PARAM_RL, args.model, MODEL_K, MODEL_THETA, MODEL_PCR, THREADS)
        align_fastq = 'bwa mem %s %s.fastq | samtools view -bS - > %s.bam'%(REFFA, output_prefix, output_prefix)
        index_bam = 'samtools sort %s.bam > %s.sorted.bam; samtools index %s.sorted.bam'%(output_prefix, output_prefix, output_prefix)
        remove_fastq = 'rm %s*.fastq'%(output_prefix)
        make_tags = "mkdir -p %s; makeTagDirectory %s %s.sorted.bam"%(output_prefix, output_prefix, output_prefix)
        find_peaks = "findPeaks %s -o auto"%(output_prefix)
        convert_to_bed = "pos2bed.pl %s/peaks.txt > %s.bed"%(output_prefix, output_prefix)

        for cmd in [chipmunk, align_fastq, index_bam, remove_fastq, make_tags, find_peaks, convert_to_bed]:
            p = Popen(cmd, shell=True)
            p.communicate()
            if p.returncode != 0: sys.exit("Error running: %s\n"%cmd)

if __name__ == "__main__":
    main()
