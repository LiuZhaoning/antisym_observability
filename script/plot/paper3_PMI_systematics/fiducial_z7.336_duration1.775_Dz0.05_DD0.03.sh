#!/bin/bash  
#PBS -N 3_PMI_SYS
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10


#z_r 7.336, Delta_z = 1.745
zeta=33.58
Tvir=113239
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='general_fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

#z_r 7.336, Delta_z = 1.775
zeta=31.83
Tvir=106666
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='general_fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

#z_r 7.336, Delta_z = 1.805
zeta=30.21
Tvir=100490
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='general_fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

#---------------------------------------------------

#z_r 7.286, Delta_z = 1.775
zeta=31.25
Tvir=108122
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='general_fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

#z_r 7.386, Delta_z = 1.775
zeta=32.42
Tvir=105208
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='general_fisher_matrix'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

