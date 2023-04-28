#!/bin/bash  
#PBS -N fiducial_sensitivity
#PBS -V 
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=5

zeta=178
Tvir=461000
Rmfp=50
SMOOTHING_SCALE=300
SMOOTHING_Pk=200
FILE_NAME='fiducial_sensitivity'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE
