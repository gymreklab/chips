#!/bin/bash

# Example:
# ./process_encode.sh https://www.encodeproject.org/files/ENCFF053NVZ/@@download/ENCFF053NVZ.bam https://www.encodeproject.org/files/ENCFF307JUJ/@@download/ENCFF307JUJ.bed.gz /storage/mgymrek/del/ K562-ELF4 Both
# ./process_encode.sh https://www.encodeproject.org/files/ENCFF673WUM/@@download/ENCFF673WUM.bam https://www.encodeproject.org/files/ENCFF016QUV/@@download/ENCFF016QUV.bed.gz /storage/mgymrek/del/ GM12878-TARDP Both

BAMURL=$1
BEDURL=$2
OUTDIR=$3
FACTOR=$4
RTYPE=$5

PICPATH=/storage/resources/source/
CHIPMUNK=chipmunk

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
wget -O ${OUTDIR}/${FACTOR}/${FACTOR}.bam ${BAMURL} || die "Could not download BAM"
wget -O ${OUTDIR}/${FACTOR}/${FACTOR}.bed.gz ${BEDURL} || die "Could not download BED"

echo "Unzipping and indexing"
gunzip -f ${OUTDIR}/${FACTOR}/${FACTOR}.bed.gz || die "Could not unzip BED"
samtools index ${OUTDIR}/${FACTOR}/${FACTOR}.bam || die "Could not index BAM file"

# Flag the BAM file
echo "MarkDuplicates"
java -jar -Xmx10G -Djava.io.tmpdir=${TMPDIR} \
    ${PICPATH}/picard.jar \
    MarkDuplicates VALIDATION_STRINGENCY=SILENT \
    I=${OUTDIR}/${FACTOR}/${FACTOR}.bam \
    O=${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam \
    M=${OUTDIR}/${FACTOR}/${FACTOR}.metrics || die "Error running mark duplicates"
samtools index ${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam || die "Error indexing dup file"


if [ "${RTYPE}" = "Paired" ] || [ "${RTYPE}" = "Both" ]; then
    $CHIPMUNK learn \
	-b ${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam \
	-p ${OUTDIR}/${FACTOR}/${FACTOR}.bed \
	-o ${OUTDIR}/${FACTOR}/${FACTOR}.paired \
	-t bed \
	--paired --output-frag-lens
fi
if [ "${RTYPE}" = "Single" ] || [ "${RTYPE}" = "Both" ]; then
    $CHIPMUNK learn \
	-b ${OUTDIR}/${FACTOR}/${FACTOR}.flagged.bam \
	-p ${OUTDIR}/${FACTOR}/${FACTOR}.bed \
	-t bed \
	-o ${OUTDIR}/${FACTOR}/${FACTOR} -c 7 --thres 100
fi
