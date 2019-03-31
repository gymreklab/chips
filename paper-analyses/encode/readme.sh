#!/bin/bash

# Get list of ENCODE accessions
./get_encode_v2.py > encode_gm12878_accs.csv

# Process histone mods on snorlax
grep H3 encode_gm12878_accs.csv > encode_gm12878_accs_HM.csv
./process_encode_snorlax.sh encode_gm12878_accs_HM.csv 5

# Process other example TFs
grep "BACH1\|BCLAF1\|CTCF\|ETV6\|IKZF1\|MEF2A||RELB\|SRF\|TARDBP\|CTCF" encode_gm12878_accs.csv > encode_gm12878_accs_TFs.csv
./process_encode_snorlax.sh encode_gm12878_accs_TFs.csv 100

# Process all on AWS - see readme_aws.sh for setup
cat encode_gm12878_accs.csv | grep -v H3 | grep -v BACH1 | grep -v BCLAF1 | grep -v CTCF | grep -v ETV6 | grep -v IKZF1 | grep -v MEF2A | grep -v RELB | grep -v SRF | grep -v TARDBP > encode_gm12878_accs_foraws.csv
./run_aws_jobs.sh encode_gm12878_accs_foraws.csv
./compile_aws_params.sh > ChIPMunk_SuppTable2.csv
./get_markdown.sh ChIPMunk_SuppTable2.csv


