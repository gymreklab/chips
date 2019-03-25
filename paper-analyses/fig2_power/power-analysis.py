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
    parser.add_argument("--optargs", help="Extra arguments to chipmunk", type=str, default="")
    parser.add_argument("--debug", help="Print commands, but don't run them", action="store_true")

    args = parser.parse_args()
    OUTDIR = args.out
    read_list = [int(item) for item in args.readnums.split(",")]

    for numreads in read_list:
        commands = []
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
                    '--thread %d %s')%(CHIPMUNK, bed_file, REFFA, output_prefix, PARAM_NC, numreads, PARAM_RL, args.model, MODEL_K, MODEL_THETA, MODEL_PCR, THREADS, args.optargs)
        commands.append(chipmunk)
        chipmunk_wce = ('%s simreads '
                        '-t wce '
                        '-f %s '
                        '-o %s.wce '
                        '--numcopies %d '
                        '--numreads %d '
                        '--readlen %d '
                        '--model %s '
                        '--gamma-frag %s,%s '
                        '--pcr_rate %f '
                        '--thread %d %s')%(CHIPMUNK, REFFA, output_prefix, PARAM_NC, numreads, PARAM_RL, args.model, MODEL_K, MODEL_THETA, MODEL_PCR, THREADS, args.optargs)
        commands.append(chipmunk_wce)
        for pref in [output_prefix, output_prefix+".wce"]:
            align_fastq = 'bwa mem -t %d %s %s.fastq | samtools view -bS - > %s.bam'%(THREADS, REFFA, pref, pref)
            index_bam = 'samtools sort --threads %d -T %s %s.bam > %s.sorted.bam; samtools index %s.sorted.bam'%(THREADS-1, pref, pref, pref, pref)
            remove_fastq = 'rm %s.fastq'%(pref)
            markdup = 'java -jar  -Xmx12G -Djava.io.tmpdir=%s $PICARD MarkDuplicates VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING I=%s.sorted.bam M=%s.metrics O=%s.flagged.bam'%(pref, pref, pref, pref)
            commands.extend([align_fastq, index_bam, remove_fastq, markdup])
        macs = 'macs2 callpeak -t %s.flagged.bam -c %s.wce.flagged.bam --name %s.%s --outdir %s'%(output_prefix, output_prefix, os.path.basename(output_prefix), numreads, os.path.dirname(output_prefix))
        commands.append(macs)

#        make_tags = "mkdir -p %s; makeTagDirectory %s %s.sorted.bam"%(output_prefix, output_prefix, output_prefix)
#        find_peaks = "findPeaks %s -o auto"%(output_prefix)
#        convert_to_bed = "pos2bed.pl %s/peaks.txt > %s.bed"%(output_prefix, output_prefix)

        for cmd in commands:
            if args.debug:
                print(cmd)
                continue
            p = Popen(cmd, shell=True)
            p.communicate()
            if p.returncode != 0: sys.exit("Error running: %s\n"%cmd)

if __name__ == "__main__":
    main()
