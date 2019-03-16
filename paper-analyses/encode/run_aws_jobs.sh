#!/bin/bash

# Run process encode on 12 example factors
SCRIPTNAME=process_encode_s3.sh
while IFS='' read -r line || [[ -n "$line" ]]; do
    bamurl=$(echo $line | cut -f 4 -d',')
    bedurl=$(echo $line | cut -f 5 -d',')
    ct=$(echo $line | cut -f 1 -d',')
    f=$(echo $line | cut -f 2 -d',')
    bamacc=$(echo $bamurl | cut -d'/' -f 5)
    bedacc=$(echo $bedurl | cut -d'/' -f 5)
    factor=${ct}_${f}_${bamacc}_${bedacc}
    rtype=$(echo $line | cut -f 3 -d',')
    cmd="aws batch submit-job \
	--job-name chipmunk-encode-${factor} \
	--job-queue chipmunk-encode \
	--job-definition chipmunk-encode:1 \
	--container-overrides 'command=[\"${SCRIPTNAME}\",\"${bamurl}\",\"${bedurl}\",\"${factor}\",\"${rtype}\"],environment=[{name=\"BATCH_FILE_TYPE\",value=\"script\"},{name=\"BATCH_FILE_S3_URL\",value=\"s3://gymreklab-awsbatch/${SCRIPTNAME}\"}]'"
    echo "${cmd}"
#    sh -c "${cmd}"
#    exit 1 # exit when debugging
done < encode_datasets_K562_GM12878_clean.csv
