#!/bin/bash  
#PBS -N fiducial_sensitivity_1.6
#PBS -V 
#PBS -lselect=1:ncpus=4:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=4

#z_r = 7.336, Delta_z = 1.8
zeta=30.48
Tvir=101500
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='fiducial_sensitivity'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE
