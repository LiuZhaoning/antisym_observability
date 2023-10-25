#!/bin/bash  
#PBS -N Pk_A_group2_384
#PBS -V
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=5
Rmfp=50
SMOOTHING_SCALE=384
SMOOTHING_Pk=384

zeta=19.14
Tvir=49100
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=22.02
Tvir=43300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

zeta=25.14
Tvir=38200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $NUMCORE

