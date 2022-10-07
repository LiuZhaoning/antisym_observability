#!/bin/bash
#PBS -N fisher_Rmfp
#PBS -V
#PBS -lselect=1:ncpus=10:mem=10gb
#PBS -j oe
#PBS -k oe

zeta=25
Tvir=30000
Rmfp=50
SMOOTHING=384
NUMCORE=10

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir 49 $SMOOTHING $NUMCORE

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py $zeta $Tvir 51 $SMOOTHING $NUMCORE



