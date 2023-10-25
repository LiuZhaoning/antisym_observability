#!/bin/bash  
#PBS -N xi_A_group1
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10
Rmfp=50

#group1
zeta=20.1
Tvir=74900
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

#it has been run successfully below
#zeta=22.3
#Tvir=68800
#python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

#zeta=24.6
#Tvir=63400
#python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE
