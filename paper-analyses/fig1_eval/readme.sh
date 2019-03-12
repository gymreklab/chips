#!/bin/bash
# /storage/mlamkin/projects/encode_data_learn/datasets/learn_files/analyzed_encode_data.csv has paths

######### Sim examples #########
# Setup windows
#./make_windows.sh 

# Learn (single/paired/frags)
./learn_factor.sh

# Simulate TF example
#factor=K562+SP1
factor=K562+H3K4me3

#./get_pcr.sh ${factor} > ${OUTDIR}/${factor}.pcr.txt
./get_frags.sh ${factor}
for nc in 10 25 50 75 10 100 1000 10000
do
    ./run_factor.sh ${factor} ${nc}
done

