#!/bin/bash

source params.sh

######### Sim examples #########

# Simulate examples
OUTDIR=/storage/mgymrek/chipmunk/fig1_eval/
nc=10
./run_factor.sh GM12878_H3K4me3_ENCFF398NET_ENCFF068UAA ${nc} 36 276590 Single
./run_factor.sh GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH ${nc} 51 623524 Single
./run_factor.sh GM12878_SRF_ENCFF387RFR_ENCFF500GHH ${nc} 101 1633899 Paired
./run_factor.sh GM12878_TARDBP_ENCFF673WUM_ENCFF016QUV ${nc} 101 652216 Paired
./run_factor.sh GM12878_BACH1_ENCFF518TTP_ENCFF866OLZ ${nc} 101 665975 Paired

#for factor in GM12878_H3K4me3_ENCFF398NET_ENCFF068UAA GM12878_H3K27ac_ENCFF097SQI_ENCFF465WTH GM12878_SRF_ENCFF387RFR_ENCFF500GHH GM12878_TARDBP_ENCFF673WUM_ENCFF016QUV
#do
#    # Convert ENCODE to tdf
#    igvtools count ${ENCDIR}/${factor}/${factor}.flagged.bam ${OUTDIR}/${factor}/${factor}.encode.tdf ${REFFA}
#done

# Evaluate cell number
#./eval_cell_number.sh
