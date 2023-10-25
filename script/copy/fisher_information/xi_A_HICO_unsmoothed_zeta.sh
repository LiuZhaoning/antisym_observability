#!/bin/bash
#PBS -N fisher_zeta
#PBS -V
#PBS -lselect=1:ncpus=10:mem=10gb
#PBS -j oe
#PBS -k oe

zeta=25
Tvir=30000
Rmfp=50
SMOOTHING=384
NUMCORE=10

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py 24 $Tvir $Rmfp $SMOOTHING $NUMCORE

python3 /home/liuzhaoning/antisym_observability/xi_A_HICO_unsmoothed.py 26 $Tvir $Rmfp $SMOOTHING $NUMCORE



