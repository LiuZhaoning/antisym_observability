#!/bin/bash  
#PBS -N xi_A_black
#PBS -V 
#PBS -lselect=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -k oe

zeta=32.5
Tvir=50000
Rmfp=50
NUMCORE=20

python3 /home/liuzhaoning/antisym_observability/antisym_unsmoothed_computation.py $zeta $Tvir $Rmfp $NUMCORE