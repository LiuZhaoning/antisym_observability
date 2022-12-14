#!/bin/bash  
#PBS -N Pk_A_blue
#PBS -V 
#PBS -lselect=1:ncpus=5:mem=5gb
#PBS -j oe
#PBS -k oe

zeta=37.5
Tvir=41000
Rmfp=50
SMOOTHING=384
NUMCORE=5

python3 /home/liuzhaoning/antisym_observability/Pk_A_model_examination.py $zeta $Tvir $Rmfp $SMOOTHING  $NUMCORE
