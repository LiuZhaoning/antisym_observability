#!/bin/bash  
#PBS -N Pk_A_cyan
#PBS -V 
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

zeta=32.5
Tvir=66700
Rmfp=50
SMOOTHING=384
NUMCORE=5

python3 /home/liuzhaoning/antisym_observability/Pk_A_model_examination.py $zeta $Tvir $Rmfp $SMOOTHING $NUMCORE
