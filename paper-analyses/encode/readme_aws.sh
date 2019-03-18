#!/bin/bash

# Note: r4.large has 2 CPU, 16GB memory, 8GB storage on /

# Set up compute queue
aws batch create-compute-environment \
    --compute-environment-name r4-large-12GB \
    --type MANAGED \
    --state ENABLED \
    --compute-resources file://r4-large.json \
    --service-role arn:aws:iam::369425333806:role/service-role/AWSBatchServiceRole

aws batch create-job-queue \
    --job-queue-name chipmunk-encode \
    --state ENABLED \
    --priority 100 \
    --compute-environment-order order=1,computeEnvironment=r4-large-12GB

aws batch register-job-definition \
    --job-definition-name chipmunk-encode \
    --type container \
    --container-properties file://chipmunk-container-properties.json

# Set up docker and AWS scripts
docker build -t gymreklab/chipmunk-run .
docker push gymreklab/chipmunk-run
aws s3 cp process_encode.sh s3://gymreklab-awsbatch/process_encode.sh
aws s3 cp process_encode_s3.sh s3://gymreklab-awsbatch/process_encode_s3.sh

# Test job on Docker
AWS_ACCESS_KEY_ID=$(cat ~/.aws/credentials  | grep id | cut -f 2 -d '=' | head -n 1 | cut -f 2 -d' ')
AWS_SECRET_ACCESS_KEY=$(cat ~/.aws/credentials  | grep secret | cut -f 2 -d '=' | head -n 1 | cut -f 2 -d' ')
BAMURL=https://www.encodeproject.org/files/ENCFF489JLD/@@download/ENCFF489JLD.bam
BEDURL=https://www.encodeproject.org/files/ENCFF489ZRF/@@download/ENCFF489ZRF.bed.gz
FACTOR=GM12878_RCOR1_ENCFF489JLD_ENCFF489ZRF
RTYPE=Single
docker run \
    -v /storage/mgymrek/del:/scratch \
    --env BATCH_FILE_TYPE="script" \
    --env BATCH_FILE_S3_URL="s3://gymreklab-awsbatch/process_encode_s3.sh" \
    --env AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
    --env AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
    -it gymreklab/chipmunk-run \
    process_encode_s3.sh ${BAMURL} ${BEDURL} ${FACTOR} ${RTYPE}

# Get full ENCODE list
./process_encode_list.sh

# Run on AWS
./run_aws_jobs.sh

# Compile all params from AWS for supp table 2
./compile_aws_params.sh > ChIPMunk_SuppTable2.csv
