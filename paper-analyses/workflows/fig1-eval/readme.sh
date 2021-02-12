############################################
# The snakemake workflow in this directory takes in an
# encode dataset and simulation method
# (chips, chipulate, or ischip) and:
# 
# 1. Learns a model
# 2. Runs chips/chipulate/ischip to simulate for various values of --numcopies
# 3. Aligns reads and flags dups
# 4. Counts reads in 1kb and 5kb windows for comparson with encode data
#
# Scatterplots vs. encode and other plots can be visualized
# in ChIPs-Figure1-Eval.ipynb
# Input datasets are read from datasets.csv
############################################

./readme-chips.sh # run sims for chips
./readme-chipulate.sh # run sims for chipulate (only single-end)
