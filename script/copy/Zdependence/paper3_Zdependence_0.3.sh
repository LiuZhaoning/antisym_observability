#!/bin/bash  
#PBS -N Zdependence
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10

# z_r = 7.036, Delta_z = 1.65
zeta=35.82
Tvir=149000
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='paper3_Zdependence'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

# z_r = 7.036, Delta_z = 1.8
zeta=27.33
Tvir=110300
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='paper3_Zdependence'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

# z_r = 7.036, Delta_z = 1.95
zeta=21.59
Tvir=82410
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='paper3_Zdependence'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE


# z_r = 7.636, Delta_z = 1.65
zeta=45.1
Tvir=126600
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='paper3_Zdependence'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

# z_r = 7.636, Delta_z = 1.8
zeta=33.96
Tvir=93530
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='paper3_Zdependence'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

# z_r = 7.636, Delta_z = 1.95
zeta=26.52
Tvir=69780
Rmfp=50
SMOOTHING_SCALE=300
FILE_NAME='paper3_Zdependence'

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $FILE_NAME $NUMCORE

SMOOTHING_Pk=200
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=250
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE

SMOOTHING_Pk=300
python3 /home/liuzhaoning/antisym_observability/Pk_A_group_class.py $zeta $Tvir $Rmfp $SMOOTHING_SCALE $SMOOTHING_Pk $FILE_NAME $NUMCORE
