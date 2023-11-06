#!/bin/bash  
#PBS -N 3_PMI
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10


# z_r = 7.336, Delta_z = 1.66
zeta=39.42
Tvir=134500
Rmfp=50
SMOOTHING_SCALE=384
FILE_NAME='general_fisher_matrix'

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE



# z_r = 7.336, Delta_z = 1.88
zeta=26.69
Tvir=86760
Rmfp=50
SMOOTHING_SCALE=384
FILE_NAME='general_fisher_matrix'

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE
