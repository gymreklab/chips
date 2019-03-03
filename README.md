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
  -p <peaks.bed> \
  -t <homer|bed>
  -o <outprefix>
```

```
chipmunk simreads \
  -p <peaks.bed \
  -f <ref.fa>
  -t <homer|bed> \
  -o <outprefix>
```

<a name="detailed"></a>
## Detailed usage

### `chipmunk learn`

### `chipmunk simreads`

<a name="formats"></a>
## File formats