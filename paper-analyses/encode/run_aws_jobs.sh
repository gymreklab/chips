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
    # First check if output file exists
    if [ "${rtype}" = "Paired" ]; then
	outfile=${factor}.paired.json
	echo "Skipping paired end for now"
	continue # process paired later... most too big and fail on AWS
    elif [ "${rtype}" = "Single" ]; then
	outfile=${factor}.json
    fi
    aws s3 ls s3://chipmunk-encode-models/${outfile} > /dev/null
    if [[ $? -eq 0 ]]; then
	echo "Found file... skipping ${factor}"
	continue
    fi

    cmd="aws batch submit-job \
	--job-name chipmunk-encode-${factor} \
	--job-queue chipmunk-encode \
	--job-definition chipmunk-encode:2 \
        --timeout 'attemptDurationSeconds=3600' \
	--container-overrides 'command=[\"${SCRIPTNAME}\",\"${bamurl}\",\"${bedurl}\",\"${factor}\",\"${rtype}\"],environment=[{name=\"BATCH_FILE_TYPE\",value=\"script\"},{name=\"BATCH_FILE_S3_URL\",value=\"s3://gymreklab-awsbatch/${SCRIPTNAME}\"}]'"
    sh -c "${cmd}"
done < encode_datasets_K562_GM12878_clean.csv
