#!/bin/bash  
#PBS -N xi_A_group1
#PBS -V 
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=5

zeta=15.85
Tvir=37100
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

zeta=18.07
Tvir=32500
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

zeta=20.50
Tvir=28600
Rmfp=50
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE


