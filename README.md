# ChIPmunk

ChIPmunk is a tool for simulating ChIP-sequencing experiments.

For questions on installation or usage, please open an issue, submit a pull request, or contact An Zheng (anz023@eng.ucsd.edu).

[Download](#download) | [Install](#install) | [Basic Usage](#usage) | [Detailed usage](#detailed) | [File formats](#formats)

<a name="download"></a>
## Download

The latest ChIPmunk release is available on the [releases page](https://github.com/gymreklab/ChIPmunk/releases).

<a name="install"></a>
## Basic Install

ChIPmunk requires the third party package [htslib](http://www.htslib.org/).

If you are installing from the tarball, type the following commands.

```
tar -xzvf ChIPmunk-X.X.tar.gz
cd ChIPmunk-X.X
./configure [--prefix=$PREFIX]
make
make install
``` 

If you do not have root access, you can set `--prefix=$HOME`, which will install `chipmunk` to `~/bin/chipmunk`.
If you get a pkg-config error, you may need to set PKG_CONFIG_PATH to a directory where it kind find `htslib.pc`.

Typing `chipmunk --help` should show a help message if ChIPmunk was successfully installed.

## Compiling from git source

To compile from git source, first make sure htslib is installed. Then run:
```
git clone https://github.com/gymreklab/ChIPmunk
cd ChIPmunk/
./reconf
./configure
make
make install
```

<a name="usage"></a>
## Basic usage

ChIPmunk is a single command line tool that contains several modules. To see available modules type:

```
chipmunk
```

The following modules are available:

* `learn`: Learn key parameters from existing ChIP-seq datasets. Learn takes in alinged reads (BAM) and peaks (BED) and outputs a model parameters file.
* `simreads`: Simulate ChIP-seq reads based on model and experimental parameters. Simulate takes in peaks (BED), model parameters (either user-specified, or learned from an existing dataset using `learn`) and outputs raw reads (FASTQ).

Basic usage is shown below. See [detailed usage](#detailed) below for more info.

```
chipmunk learn \
  -b <reads.bam> \
  -p <peaks> \
  -t <homer|bed>
  -o <outprefix>
```

```
chipmunk simreads \
  -p <peaks> \
  -f <ref.fa>
  -t <homer|bed> \
  -o <outprefix>
```

<a name="detailed"></a>
## Detailed usage

### chipmunk learn

Required parameters:
* `-b <file.bam>`: BAM file containing aligned reads. Must be sorted and indexed. Paired end or single end data are supported.
* `-p <peaks>`: file containing peaks. 
* `-t <homer|bed>`: Specify the format of the peaks file. Options are "bed" or "homer".
* `-o <outprefix>`: Prefix to name output files. Outputs file `<outprefix>.json` with learned model parameters.

### chipmunk simreads

Required parameters:
* `-p <peaks>`: file containing peaks. 
* `-t <homer|bed>`: Specify the format of the peaks file. Options are "bed" or "homer".
* `-f <ref.fa>`: Reference genome fasta file. Must be indexed (e.g. `samtools faidx <ref.fa>`)
* `-o <outprefix>`: Prefix to name output files. Outputs `<outprefix>.fastq` for single-end data or `<outprefix>_1.fastq` and `<outprefix>_2.fastq` for paired-end data.

Experiment parameters:
* `--numcopies <int>`: Number of copies of the genome to simulate (Default: 100)
* `--numreads <int>`: Number of reads (or read pairs) to simulate (Default: 1000000)
* `--readlen <int>`: Read length to generate (Default: 36bp)
* --paired`: Simulated paired-end reads (by default single-end reads are generated).

Model parameters: (either user-specified or learned from `chipmunk learn`:
* `--model <str>`: JSON file with model parameters (e.g. from running learn. Setting parameters with other options overrides anything in the JSON file.
* `--gamma-frag <float>,<float>`: Parameters for fragment length distribution (k, theta for Gamma distribution). Default: 15.67,15,49
* `--spot <float>`: SPOT score (fraction of reads in peaks). Default: 0.18
* `--frac <float>`: Fraction of the genome that is bound. Default: 0.03
* `--pcr_rate <float>`: The geometric step size paramters for simulating PCR. Default: 1.0.

Peark scoring:
* `-b <reads.bam>`: Use a provided BAM file to obtain scores for each peak. No BAM is required. If a BAM is not given, scores in the peak files are used.
* `-c <int>`: The index of the BED file column used to score each peak (index starting from 1). Required if not using `-b`.

Other options:
* `--region <str>`: Only simulate reads from this region chrom:start-end. By default, simulate genome-wide.
* `--binsize <int>`: Consider bins of this size when simulating. Default: 100000.
* `--thread <int>`: Number of threads to use. Default: 1.
* `--sequencer <str>`: Sequencing error mode. If not set, use `--sub`,`--ins`, and `--del`. Specify `--sequencer HiSeq` to set `--sub 2.65e-3 --del 2.43e-4 --ins 1.83e-4`.
* `--sub <float>`: Substitution error rate. Default: 0.
* `--ins <float>`: Insertion error rate. Default: 0.
* `--del <float>`: Deletion error rate. Default: 0.

<a name="formats"></a>
## Formats

### Peak files

Peak files may be either in Homer or Bed format. For all modules, the option `-t` should specify either "homer" or "bed" appropriately.

With `-t bed`, your peak file must be tab-delimited with no header line. The first three columns are chromosome, start, and end. For `chipmunk simreads` if you don't supply a BAM file you'll also need to specify which column contains the peak score using `-c <colnum>`. For example if your file just has four columns, chrom, start, end, and score, set `-c 4`. If your peaks don't have scores you can modify your peak file to set all peaks to have score 1.

With `-t homer`, your peak file should be in the format output by the [Homer peak caller](http://homer.ucsd.edu/homer/ngs/peaks.html). For homer files, set `-c 6` since the peak score intensity is in column 6.

### Model files

Model files are in JSON syntax, and follow the example below:

```
{
    "frag": {
        "k": 14.176023483276367,
        "theta": 16.013242721557617
    },
    "pcr_rate": 0.825218915939331,
    "pulldown": {
        "f": 0.043268248438835144,
        "s": 0.1925637274980545
    }
}
```

`chipmunk learn` outputs a JSON model file. `chipmunk simreads` can take in a model file with all or some of these parameters specified. Model parameters set on the command line override those set in the JSON model file.