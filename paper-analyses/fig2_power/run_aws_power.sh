#!/bin/bash

aws s3 cp chipmunk_power_s3.sh s3://gymreklab-awsbatch/chipmunk_power_s3.sh
SCRIPTNAME=chipmunk_power_s3.sh

factors="GM12878_H3K27me3_ENCFF014HHB_ENCFF533IKU GM12878_H3K36me3_ENCFF191SDM_ENCFF695NNX GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ GM12878_BACH1_ENCFF518TTP_ENCFF866OLZ GM12878_CTCF_ENCFF584BRF_ENCFF559IXF"
readnums="1000000 5000000 10000000 25000000 50000000 100000000"

for factor in $factors
do
    for numreads in $readnums
    do
	cmd="aws batch submit-job \
    --job-name chipmunk-power-${factor}-${numreads} \
    --job-queue chipmunk-power \
    --job-definition chipmunk-power:1 \
    --timeout 'attemptDurationSeconds=100000' \
    --container-overrides 'command=[\"${SCRIPTNAME}\",\"${factor}\",\"${numreads}\"],environment=[{name=\"BATCH_FILE_TYPE\",value=\"script\"},{name=\"BATCH_FILE_S3_URL\",value=\"s3://gymreklab-awsbatch/${SCRIPTNAME}\"}]'"
	sh -c "${cmd}"
    done
done
