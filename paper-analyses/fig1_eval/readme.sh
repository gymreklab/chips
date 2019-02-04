#!/bin/bash

######### Sim examples #########
# Setup windows
./make_windows.sh 

# Simulate TF example
factor=K562+SP1
#./get_frags.sh ${factor}
for nc in 25 50 75 #10 100 1000
do
    ./run_factor.sh ${factor} ${nc}
done

