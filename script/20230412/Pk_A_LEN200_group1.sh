#!/bin/bash  
#PBS -N Pk_A_LEN200_group1
#PBS -V 
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=5
Rmfp=50
SMOOTHING_SCALE=384
SMOOTHING_Pk=200

zeta=20.1
Tvir=74900
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=22.3
Tvir=68800
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=24.6
Tvir=63400
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE


