# Set up compute

```
# Build docker
docker build -t gymreklab/chipmunk-power-run .
docker push gymreklab/chipmunk-power-run

# Make queue
aws batch create-job-queue \
    --job-queue-name chipmunk-power \
    --state ENABLED \
    --priority 100 \
    --compute-environment-order order=1,computeEnvironment=c4-2xlarge-100GB

# Use queue chipmunk-power, job definition chipmunk-power
aws batch register-job-definition \
    --job-definition-name chipmunk-power \
    --type container \
    --container-properties file://chipmunk-power-container-properties.json


```

# Upload ENCODE data for factors we want
```
factors="GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH GM12878_H3K4me3_ENCFF398NET_ENCFF068UAA GM12878_H3K27me3_ENCFF014HHB_ENCFF533IKU GM12878_H3K36me3_ENCFF191SDM_ENCFF695NNX GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ GM12878_BACH1_ENCFF518TTP_ENCFF866OLZ GM12878_SRF_ENCFF387RFR_ENCFF500GHH"
for factor in $factors
do
    aws s3 cp /storage/mgymrek/chipmunk/encode/${factor}/${factor}.bed s3://gymreklab/chipmunk-power/${factor}/${factor}.bed
    aws s3 cp /storage/mgymrek/chipmunk/encode/${factor}/${factor}.json s3://gymreklab/chipmunk-power/${factor}/${factor}.json
done
```


# Test job in docker
aws s3 cp chipmunk_power_s3.sh s3://gymreklab-awsbatch/chipmunk_power_s3.sh
AWS_ACCESS_KEY_ID=$(cat ~/.aws/credentials  | grep id | cut -f 2 -d '=' | head -n 1 | cut -f 2 -d' ')
AWS_SECRET_ACCESS_KEY=$(cat ~/.aws/credentials  | grep secret | cut -f 2 -d '=' | head -n 1 | cut -f 2 -d' ')
factor=GM12878_H3K4me3_ENCFF398NET_ENCFF068UAA
numreads=1000000
optargs="--region chr19:1-59128983"
docker run \
    -v /storage/mgymrek/del:/scratch \
    --env BATCH_FILE_TYPE="script" \
    --env BATCH_FILE_S3_URL="s3://gymreklab-awsbatch/chipmunk_power_s3.sh" \
    --env AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
    --env AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
    -it gymreklab/chipmunk-power-run \
    chipmunk_power_s3.sh ${factor} ${numreads} "${optargs}"

# Test job on AWS
```
aws s3 cp chipmunk_power_s3.sh s3://gymreklab-awsbatch/chipmunk_power_s3.sh
factor=GM12878_H3K27me3_ENCFF014HHB_ENCFF533IKU
numreads=1000000
SCRIPTNAME=chipmunk_power_s3.sh
cmd="aws batch submit-job \
    --job-name chipmunk-power-${factor}-${numreads} \
    --job-queue chipmunk-power \
    --job-definition chipmunk-power:1 \
    --timeout 'attemptDurationSeconds=86400' \
    --container-overrides 'command=[\"${SCRIPTNAME}\",\"${factor}\",\"${numreads}\"],environment=[{name=\"BATCH_FILE_TYPE\",value=\"script\"},{name=\"BATCH_FILE_S3_URL\",value=\"s3://gymreklab-awsbatch/${SCRIPTNAME}\"}]'"
echo $cmd
```

# Run on AWS
./run_aws_power.sh

# Get AWS peak results
aws s3 ls s3://gymreklab/chipmunk-power/ |grep Peak | grep "03-29" | awk '{print $NF}' | \
    xargs -n1 -P1 -I% sh -c "aws s3 cp s3://gymreklab/chipmunk-power/% /storage/mgymrek/chipmunk/fig2_power/froms3/"

# Get Power with AWS peaks
```
factors="GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH GM12878_H3K4me3_ENCFF398NET_ENCFF068UAA GM12878_H3K27me3_ENCFF014HHB_ENCFF533IKU GM12878_H3K36me3_ENCFF191SDM_ENCFF695NNX GM12878_H3K4me1_ENCFF252ZII_ENCFF966LMJ GM12878_BACH1_ENCFF518TTP_ENCFF866OLZ GM12878_CTCF_ENCFF584BRF_ENCFF559IXF"

for factor in $factors
do
  ./summ_power_aws.sh ${factor}
  ./summ_power_broad.sh ${factor}
done
```
