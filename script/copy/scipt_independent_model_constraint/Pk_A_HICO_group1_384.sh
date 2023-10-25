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

zeta=15.85
Tvir=37100
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=18.07
Tvir=32500
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=20.50
Tvir=28600
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE


