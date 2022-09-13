#!/bin/bash  
#PBS -N ACC_z28.4_T66000
#PBS -V 
#PBS -lselect=1:ncpus=10:mem=10gb
#PBS -j oe
#PBS -k oe

zeta=28.4
Tvir=66000
NUMCORE=10

if python3 /home/liuzhaoning/antisym/antisym_parameters/BMF_calculation/BMF_calculation.py $zeta $Tvir $NUMCORE
then
    echo "BMF calculation finishes"
    fi
if python3 /home/liuzhaoning/antisym/antisym_parameters/HIrho_calculation/HIrho_calculation.py $zeta $Tvir $NUMCORE
then
    echo "HIrho calculation finishes"
    fi
if python3 /home/liuzhaoning/antisym/antisym_parameters/xi_A_HICO_unsmoothed_calculation/xi_A_HICO_unsmoothed_calculation.py $zeta $Tvir $NUMCORE
then
    echo "xi_A_HICO calculation finishes"
    fi
if python3 /home/liuzhaoning/antisym/antisym_parameters/Pk_A/Pk_A_computation.py $zeta $Tvir $NUMCORE
then
    echo "xi_A calculation finishes"
    fi
