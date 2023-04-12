#!/bin/bash
#PBS -N Pk_A_group3_200
#PBS -V
#PBS -lselect=1:ncpus=10:mem=3gb
#PBS -j oe
#PBS -k oe

NUMCORE=10
Rmfp=50
SMOOTHING_SCALE=384
SMOOTHING_Pk=200

zeta=23.19
Tvir=63800
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=26.82
Tvir=56200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=30.93
Tvir=49800
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE


