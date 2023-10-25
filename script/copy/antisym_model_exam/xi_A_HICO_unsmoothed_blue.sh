#!/bin/bash  
#PBS -N xi_A_blue
#PBS -V 
#PBS -lselect=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -k oe

zeta=37.5
Tvir=41000
Rmfp=50
SMOOTHING=384
NUMCORE=20

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir $Rmfp $SMOOTHING  $NUMCORE
