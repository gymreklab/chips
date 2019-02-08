# asimon-chip-sim
Simulation tool for ChIP- and other -seq experiments



## Main functions 
### learn
### Structure
```
asimon learn -b BAM_file.bam -p Peaks.bed -t Peak_File_type -o /path/to/output/prefix_to_file
```
### simreads
### Structure
```
asimon simreads 
    -p Peaks.bed \
    -t TYPE \
    -f genome_reference.fa \
    -o /path/to/output/prefix_to_file \
    -b BAM_file.bam \
    --numcopies # copies \
    --numreads # reads \
    --readlen length of reads \
    --gamma-frag alpha,beta \ 
    --spot spot_score --frac frac_score \
    --region chr:start-end \
    --binsize bin_length \
    --thread # threads \
    --sequencer Type of sequencer\
    --pcr_rate rate \
    --paired

```

##Simreads default parameters:
```
numcopies = 100
numreads = 1000000
readlen = 36
paired = False
gamma-frag = alpha: 15.67, beta: 15.49
spot = 0.17594
frac = 0.03713
region = "" (Whole Genome)
binsize = 100000
thread = 1
pcr rate = <Need to fill out>
```

Note: For the gamma distribution parameters alpha = k, beta = theta.

