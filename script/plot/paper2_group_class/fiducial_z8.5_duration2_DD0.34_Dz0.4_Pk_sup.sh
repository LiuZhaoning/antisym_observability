#!/bin/bash  
#PBS -N 2_group
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10

# z_r = 8.5, Delta_z = 2
zeta=32.71
Tvir=50269
Rmfp=50
SMOOTHING_SCALE=384
FILE_NAME='paper2_group_class'

SMOOTHING_Pk=96
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

# z_r = 8.9, Delta_z = 2
zeta=37.27
Tvir=45260
Rmfp=50
SMOOTHING_SCALE=384
FILE_NAME='paper2_group_class'

SMOOTHING_Pk=96
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE
