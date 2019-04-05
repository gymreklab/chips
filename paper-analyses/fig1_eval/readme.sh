#!/bin/bash

source params.sh

######### Sim examples #########
./eval_cell_number.sh 5 GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH 51 269113 Single
./eval_cell_number.sh 1 GM12878_CTCF_ENCFF406XWF_ENCFF833FTF 36 521557 Single

# Simulate examples
# ./sim_examples.sh

# Evaluate cell number
#./eval_cell_number.sh
