#!/bin/bash  
#PBS -N Pk_A_group1_384
#PBS -V 
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=5
Rmfp=50
SMOOTHING_SCALE=384
SMOOTHING_Pk=384

zeta=17.15
Tvir=27100
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=20.3
Tvir=22800
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=24
Tvir=19000
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE


