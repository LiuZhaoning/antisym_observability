#!/bin/bash  
#PBS -N xi_A_group2
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=5gb
#PBS -j oe
#PBS -k oe

NUMCORE=10
Rmfp=50

#group2
zeta=62.3
Tvir=252000
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

zeta=70.8
Tvir=232000
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE

zeta=80.5
Tvir=215000
python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp 384 $NUMCORE
