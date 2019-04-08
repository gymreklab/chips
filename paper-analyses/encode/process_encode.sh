#!/bin/bash

# Example:
# ./process_encode.sh https://www.encodeproject.org/files/ENCFF053NVZ/@@download/ENCFF053NVZ.bam https://www.encodeproject.org/files/ENCFF307JUJ/@@download/ENCFF307JUJ.bed.gz /storage/mgymrek/del/ K562-ELF4 Both
# ./process_encode.sh https://www.encodeproject.org/files/ENCFF673WUM/@@download/ENCFF673WUM.bam https://www.encodeproject.org/files/ENCFF016QUV/@@download/ENCFF016QUV.bed.gz /storage/mgymrek/del/ GM12878-TARDP Both

BAMURL=$1
BEDURL=$2
OUTDIR=$3
FACTOR=$4
RTYPE=$5
THRESH=$6


#MAXFILESIZE=4294967296 # 4GB. Skip large files so AWS doesn't run out of space. 
MAXFILESIZE=4294967296000

CHIPMUNK=chipmunk #/home/mgymrek/workspace/ChIPmunk/src/chipmunk

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# Set up
mkdir -p ${OUTDIR}/${FACTOR}/tmp || die "Could not make ${OUTDIR}/${FACTOR}"
TMPDIR=${OUTDIR}/${FACTOR}/tmp

# Download ENCODE data
echo "Downloading ENCODE data"
if [ ! -f ${OUTDIR}/${FACTOR}/${FACTOR}.bam ]; then
    curl -s -L --max-filesize ${MAXFILESIZE} -o ${OUTDIR}/${FACTOR}/${FACTOR}.bam ${BAMURL} || die "Could not download BAM"
    samtools index ${OUTDIR}/${FACTOR}/${FACTOR}.bam || die "Could not index BAM file"
fi

if [ ! -f ${OUTDIR}/${FACTOR}/${FACTOR}.bed ]; then
    curl -s -L --max-filesize ${MAXFILESIZE} -o ${OUTDIR}/${FACTOR}/${FACTOR}.bed.gz ${BEDURL} || die "Could not download BED"
    gunzip -f ${OUTDIR}/${FACTOR}/${FACTOR}.bed.gz || die "Could not unzip BED"
fi


# Flag the BAM file
# Need $PICARD env var set
echo "MarkDuplicates"
java -jar -Xmx12G -Djava.io.tmpdir=${TMPDIR} $PICARD \
    MarkDuplicates VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING \
    I=${OUTDIR}/${FACTOR}/${FACTOR}.bam \
    O=${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam \
    M=${OUTDIR}/${FACTOR}/${FACTOR}.metrics || die "Error running mark duplicates"
samtools index ${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam || die "Error indexing dup file"

if [ "${RTYPE}" = "Paired" ] || [ "${RTYPE}" = "Both" ]; then
    $CHIPMUNK learn \
	-b ${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam \
	-p ${OUTDIR}/${FACTOR}/${FACTOR}.bed \
	-o ${OUTDIR}/${FACTOR}/${FACTOR}.paired-1.9 \
	-t bed --scale-outliers \
	--paired --output-frag-lens
fi
if [ "${RTYPE}" = "Single" ] || [ "${RTYPE}" = "Both" ]; then
    $CHIPMUNK learn \
	-b ${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam \
	-p ${OUTDIR}/${FACTOR}/${FACTOR}.bed \
	-t bed --scale-outliers \
	-o ${OUTDIR}/${FACTOR}/${FACTOR}-1.9 -c 7 --thres $THRESH
fi
