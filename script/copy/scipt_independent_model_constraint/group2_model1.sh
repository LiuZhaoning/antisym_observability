#!/bin/bash  
#PBS -N group2_model1
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=10gb
#PBS -j oe
#PBS -k oe

NUMCORE=10

zeta=19.14
Tvir=49100
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

SMOOTHING_SCALE=384
SMOOTHING_Pk=384
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

