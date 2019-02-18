#asimon-chip-sim
Simulation tool for ChIP- and other -seq experiments


##Usage
asimon reads [options]

###Required arguments:
1. -p <peaks>:
    The path to your peak file.
    Currently, ChIPMunk support BED format, and homer format. You need to specify the format with "-t".
2. -t <str>:
    The file format of your peak file.
    You can choose from "bed" and "homer".
3. -f <ref.fa>
    The path to your reference genome FASTA file
4. -f <outprefix>
    Prefix for output files

###Experiment parameters:
5. --numcopies <int>
    Number of copies of the genome to simulate
6. --numreads <int>
    Number of reads (or read pairs) to simulate
7. --readlen <int>
    Read length to generate
8. --paired
    Simulate paired-end reads

###Model parameters
9. --gamma-frag <float>,<float>
    Parameters for fragment length distribution (alpha, beta).
