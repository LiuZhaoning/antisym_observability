#!/bin/bash  
#PBS -N xi_A_group1
#PBS -V 
#PBS -lselect=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -k oe

NUMCORE=20

zeta=17.15
Tvir=27100
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

zeta=20.3
Tvir=22800
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

zeta=24
Tvir=19000
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE


