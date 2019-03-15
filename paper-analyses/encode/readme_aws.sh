#!/bin/bash

# Note: r5.large has 2 CPU, 16GB memory, 8GB storage on /
# TODO: make chipmunk-1.7-run docker after fix seg fault

# Set up compute queue
aws batch create-compute-environment \
    --compute-environment-name r5-large \
    --type MANAGED \
    --state ENABLED \
    --compute-resources file://r5-large.json \
    --service-role arn:aws:iam::369425333806:role/service-role/AWSBatchServiceRole

aws batch create-job-queue \
    --job-queue-name chipmunk-encode \
    --state ENABLED \
    --priority 100 \
    --compute-environment-order order=1,computeEnvironment=r5-large

aws batch register-job-definition \
    --job-definition-name chipmunk-encode \
    --type container \
    --container-properties file://chipmunk-container-properties.json

./run_aws_jobs.sh
