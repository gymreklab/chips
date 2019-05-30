# Tulip

Tulip is a tool for simulating ChIP-sequencing experiments.

For questions on installation or usage, please open an issue, submit a pull request, or contact An Zheng (anz023@eng.ucsd.edu).

[Download](#download) | [Install](#install) | [From Docker](#docker) | [Basic Usage](#usage) | [Detailed usage](#detailed) | [File formats](#formats) | [FAQ](#faq)

<a name="download"></a>
## Download

The latest Tulip release is available on the [releases page](https://github.com/gymreklab/tulip/releases).

<a name="install"></a>
## Installation
### Install from Tarball
Note: Tulip uses a helper tool "pkg-config" in compilation, and requires a third party package [htslib](http://www.htslib.org/).

You can use the following commands to install Tulip from the tarball file:
```
tar -xzvf tulip-X.X.tar.gz
cd tulip-X.X
./configure
make
make install
``` 

If you do not have root access, you can install Tulip in a local directory with the commands:
```
./configure --prefix=$YOUR_PATH
make
make install
```
which installs `Tulip` to `$YOUR_PATH/bin/tulip`.

Typing `tulip --help` should show a help message if Tulip is successfully installed.

If you get a pkg-config error, you may need to set the environment variable PKG_CONFIG_PATH to the directory containing `htslib.pc`. If you do not have htslib installed in your machine, you need to install htslib first. You can download the package from https://github.com/samtools/htslib/releases/, and install htslib with the following commands:
```
tar -xjvf htslib-1.8.tar.bz2
cd htslib-1.8
./configure --disable-lzma --disable-bz2 --prefix=$YOUR_PATH
make
make install
```
Once you succussfully install the htslib, you should set the environment variables `PKG_CONFIG_PATH` to the directory containing `htslib.pc` and `LD_LIBRARY_PATH` to the directory containing the `libhts.so*` files. e.g.:
```
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:$YOUR_PATH/lib/pkgconfig
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$YOUR_PATH/lib
```

### Install from git source

To compile from git source, first make sure htslib is installed. Then run:
```
git clone https://github.com/gymreklab/tulip
cd tulip/
./reconf
./configure
make
make install
```

If installing from source on OSX, you will need to have autotools and related packages installed (ideally using Homebrew):
```
brew install autoconf automake autoconf-archive libtool pkg-config
```

<a name="docker"></a>
## From docker

Tulip is also packaged in a docker container available on the [Gymrek Lab docker hub](https://hub.docker.com/u/gymreklab). The latest release can be found at `gymreklab/tulip-1.10.1` built from the `Dockerfile` in this repo. Example usage below:

```
# Pull the image
docker pull gymreklab/tulip-1.10.1
# Test that it works
docker run gymreklab/tulip-1.10.1 tulip learn --help
```

<a name="usage"></a>
## Basic usage

Tulip is a single command line tool that contains several modules. To see available modules type:

```
tulip
```

The following modules are available:

* `learn`: Learn key parameters from existing ChIP-seq datasets. Learn takes in alinged reads (BAM) and peaks (BED) and outputs a model parameters file.
* `simreads`: Simulate ChIP-seq reads based on model and experimental parameters. Simulate takes in peaks (BED), model parameters (either user-specified, or learned from an existing dataset using `learn`) and outputs raw reads (FASTQ).

Basic usage is shown below. See [detailed usage](#detailed) below for more info.

```
tulip learn \
  -b <reads.bam> \
  -p <peaks> \
  -t <homer|bed>
  -o <outprefix>
```

```
tulip simreads \
  -p <peaks> \
  -f <ref.fa>
  -t <homer|bed|wce> \
  -o <outprefix>
```

<a name="detailed"></a>
## Detailed usage

### tulip learn

Required parameters:
* `-b <file.bam>`: BAM file containing aligned reads. Should be sorted and indexed. To accurately estimate PCR duplicate rate, duplicates must be flagged e.g. using Picard. Both paired-end or single-end data are supported.
* `-p <peaks>`: file containing peaks. 
* `-t <homer|bed>`: Specify the format of the peaks file. Options are "bed" or "homer".
* `-c <int>`: The index of the BED or homer peak file column used to score each peak (index starting from 1)
* `-o <outprefix>`: Prefix to name output files. Outputs file `<outprefix>.json` with learned model parameters.

Optional parameters for BAM parsing:
* `--paired`: Data is paired
* `--thres <float>`: For estimating fragment length distribution from single end data, only consider peaks with scores above this threshold.

Other optional parameters:
* `--scale-outliers`: Set all peaks with scores $>$2x median score to have binding prob 1. Recommended with real data.
* `--noscale`: Don't scale peak scores. Treat given scores as binding probabilities.
* `--est <int>`: Estimated fragment length. Used as a rough guess to guide inference of fragment length distribution from single end data.
* `-r <float>`: Ignore peaks with top r% of peak scores.

### tulip simreads

Required parameters:
* `-p <peaks>`: file containing peaks. 
* `-t <homer|bed|wce>`: Specify the format of the peaks file. Options are "bed" or "homer" when loading peaks. Specify `-t wce` and no peaks input file to simulate whole cell extract control data.
* `-f <ref.fa>`: Reference genome fasta file. Must be indexed (e.g. `samtools faidx <ref.fa>`)
* `-o <outprefix>`: Prefix to name output files. Outputs `<outprefix>.fastq` for single-end data or `<outprefix>_1.fastq` and `<outprefix>_2.fastq` for paired-end data.

Experiment parameters:
* `--numcopies <int>`: Number of simulation rounds (copies of the reference genome) to simulate (Default: 100). Note, this is not directly comparable to the number of input cells.
* `--numreads <int>`: Number of reads (or read pairs) to simulate (Default: 1000000)
* `--readlen <int>`: Read length to generate (Default: 36bp)
* `--paired`: Simulated paired-end reads (by default single-end reads are generated).

Model parameters: (either user-specified or learned from `tulip learn`:
* `--model <str>`: JSON file with model parameters (e.g. from running learn. Setting parameters with other options overrides anything in the JSON file.
* `--gamma-frag <float>,<float>`: Parameters for fragment length distribution (k, theta for Gamma distribution). Default: 15.67,15,49
* `--spot <float>`: SPOT score (fraction of reads in peaks). Default: 0.18
* `--frac <float>`: Fraction of the genome that is bound. Default: 0.03
* `--pcr_rate <float>`: The geometric step size paramters for simulating PCR. Default: 1.0.
* `--recomputeF`: Recompute `--frac` param based on input peaks. Recommended especially when using model parameters that were not learned on real data.

Peak scoring:
* `-b <reads.bam>`: Use a provided BAM file to obtain scores for each peak (optional). If a BAM is not given, scores in the peak files are used.
* `-c <int>`: The index of the BED or homer peak file column used to score each peak (index starting from 1). Required if not using `-b`.
* `--scale-outliers`: Set all peaks with scores >2x median score to have binding prob 1. Recommended with real data.
* `--noscale`: Don't scale peak scores. Treat given scores as binding probabilities.

Other options:
* `--seed <unsigned>`: The random seed used for initiating randomization opertions. By default or 0, use wall-clock time.
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

Peak files may be either in BED or HOMER peak format. For all modules, the option `-t` should specify either "bed" or "homer" appropriately.

With `-t bed`, your peak file must be tab-delimited with no header line. The first three columns are chromosome, start, and end. For `tulip simreads` if you don't supply a BAM file you'll also need to specify which column contains the peak score using `-c <colnum>`. For example if your file just has four columns, chrom, start, end, and score, set `-c 4`. If your peaks don't have scores you can modify your peak file to set all peaks to have score 1.

With `-t homer`, your peak file should be in the format output by the [HOMER peak caller](http://homer.ucsd.edu/homer/ngs/peaks.html). For HOMER files, set `-c 6` since the peak score intensity is in column 6.

### Model files

Model files are in JSON syntax, and follow the example below. Hundreds of model files trained on ENCODE ChIP-seq datasets for GM12878 are available for download on the [Tulip wiki](https://github.com/gymreklab/tulip/wiki/Tulip-model-files-for-GM12878-ENCODE-Datasets).

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

`tulip learn` outputs a JSON model file. `tulip simreads` can take in a model file with all or some of these parameters specified. Model parameters set on the command line override those set in the JSON model file. 

<a name="faq"></a>
## FAQ

**Q**: What should I set the number of genome copies (`--numcopies`) parameter to for `simreads`?<br>
**A**: This gives the number of simulation rounds to perform. This number is not directly comparable to the actual number of cells used in an experiment since we do not currently model pulldown inefficiency. We have found that for histone modifications performance starts to plateau after around 25 copies (`--numcopies 25`). For transcription factors we recommend setting `--numcopies 1000`. Note, run time increases linearly with the value set for this parameter.
<br><br>
**Q**: I get a `pcr_rate` output from `learn` of 1.0 (no duplicates) but I know there should be duplicates in my data!<br>
**A**: Make sure duplicates are marked, e.g. using [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).
<br><br>
**Q**: What should I do if I want to replicate my simulation experiment?<br>
**A**: Each time you run Tulip, it prints in your console the random seed being used. If you want to replicate this simulation experiment, you can simply set up "--seed" option in the simreads module with that random seed. Note: if you are using multi-threads in your work, the order of output reads in the fastq file may be different in each repeated experiment, because threads claim their jobs in an arbitrary manner during run-time. The total contents, however, stay the same in different runs.

