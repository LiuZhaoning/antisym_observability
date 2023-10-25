#!/bin/bash
#PBS -N fisher_Tvir
#PBS -V
#PBS -lselect=1:ncpus=10:mem=10gb
#PBS -j oe
#PBS -k oe

zeta=25
Tvir=30000
Rmfp=50
SMOOTHING=384
NUMCORE=10

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta 29000 $Rmfp $SMOOTHING $NUMCORE

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta 31000 $Rmfp $SMOOTHING $NUMCORE



