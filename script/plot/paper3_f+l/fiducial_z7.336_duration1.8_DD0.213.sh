#!/bin/bash  
#PBS -N f+l
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10


# z_r = 7.336, Delta_z = 1.8 - 0.213
zeta=45.73
Tvir=156300
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE



# z_r = 7.336, Delta_z = 1.8 + 0.213
zeta=21.82
Tvir=67230
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE
