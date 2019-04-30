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
cat ../encode/encode_gm12878_accs.csv | \
    grep "ENCFF385RWJ\|ENCFF398NET\|ENCFF809FFP\|ENCFF677IFV\|ENCFF191SDM\|ENCFF518TTP\|ENCFF406XW" > \
    encode_power.csv
factors=$(cat encode_power.csv  | awk -F"," '{print $1"_"$2 "/"$4"/"$5}' | cut -d'/' -f 1,8,15 | sed 's/\//_/g' | sed 's/.bam//' | sed 's/.bed.gz//')

for factor in $factors
do
    aws s3 cp /storage/mgymrek/chipmunk_round2/encode/${factor}/${factor}.bed s3://gymreklab/chipmunk-power/${factor}/${factor}.bed
    aws s3 cp /storage/mgymrek/chipmunk_round2/encode/${factor}/${factor}-1.9.json s3://gymreklab/chipmunk-power/${factor}/${factor}.json
done
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
./run_aws_power.sh $factors

# Get AWS peak results
aws s3 ls s3://gymreklab/chipmunk-power/ |grep Peak | grep "04-19" | awk '{print $NF}' | \
    xargs -n1 -P1 -I% sh -c "aws s3 cp s3://gymreklab/chipmunk-power/% /storage/mgymrek/chipmunk_round2/fig2_power/froms3/"

# Get Power with AWS peaks
```
for factor in $factors
do
  ./summ_power_aws.sh ${factor}
  ./summ_power_broad_aws.sh ${factor}
done
```
