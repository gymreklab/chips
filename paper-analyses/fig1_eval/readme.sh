#!/bin/bash

######### Sim examples #########
# Setup windows on chr19
#./make_windows.sh 

# Simulate examples
OUTDIR=/storage/mgymrek/chipmunk/fig1_eval/
for factor in GM12878_TARDBP_ENCFF673WUM_ENCFF016QUV
do
    ./run_factor.sh ${factor} ${nc}
done

