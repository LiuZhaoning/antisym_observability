import numpy as np
import math
import matplotlib 
from matplotlib import pyplot as plt
from scipy import integrate
from scipy import special
from scipy.interpolate import interp1d, interp2d
from scipy.misc import derivative
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import time
import sys
import mpmath
import os
import antisym_func
from multiprocessing import Pool
import Pk_A_function

if __name__ == '__main__':
    zeta = float(sys.argv[1])
    T_vir = float(sys.argv[2])
    R_mfp = float(sys.argv[3])
    NUM_CORE = int(sys.argv[4])
    BOX_LEN = 300 
    eta = 0.5
    M_max = antisym_func.RtoM(R_mfp)
    if (T_vir < 9.99999e3):
        mu = 1.22
    else:
        mu = 0.6
        
    #calculate some of the zeta_z points to do interpolation
    NUM = 100
    z_zeta_interp_array = np.linspace(5.5, 13.5, NUM)
    zeta_z_interp_array = []
    for i in range(NUM):
        zeta_z_interp_array.append(antisym_func.zeta_z(z_zeta_interp_array[i], zeta, T_vir, mu))
    #zeta(z) interpolation
    zeta_z_func = interp1d(z_zeta_interp_array, zeta_z_interp_array, kind = 'cubic')
    
    #calculate the history
    z_history = np.linspace(6, 12, 100); HI_history = []
    for z in z_history:
        HI_history.append(1 - antisym_func.bar_Q(z, M_max, zeta_z_func, T_vir, mu, antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)))
    [z_dxHdz_history, dxHdz_history] = antisym_func.dxH_dz_cal(z_history, HI_history)
    max_dxHdz = max(dxHdz_history)

    #load in the unsmoothed antisymmetric cross-correlation data
    xi_A_HICO_unsmoothed_map = [] #for this condition specifically
    z_grid = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g/z_grid'%(zeta, T_vir, R_mfp))
    r12_grid = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g/r12_grid'%(zeta, T_vir, R_mfp))
    for z in z_grid:
        xi_A_HICO_unsmoothed_map.append(np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g/xi_A_HICO_unsmoothed_z%4.4g.txt'%(zeta, T_vir, R_mfp, z)))
    xi_A_HICO_unsmoothed_func = interp2d(r12_grid, z_grid, xi_A_HICO_unsmoothed_map, kind = 'cubic')
        
    #set the redshift of power spectrum we are computating
    dxHdz_temp_array = np.linspace(0.225, 0.44, 25)
    z_xi_acc_smoothed_array = []; z_xi_dec_smoothed_array = [];
    dxHdz_xi_acc_smoothed_array = []; dxHdz_xi_dec_smoothed_array = []
    for dxHdz in dxHdz_temp_array:
        if (dxHdz < max_dxHdz):
            dxHdz_xi_acc_smoothed_array.append(dxHdz)
            z_xi_acc_smoothed_array.append(antisym_func.dxHdz_To_z(dxHdz, M_max, zeta_z_func, T_vir, mu, z_dxHdz_history, dxHdz_history)[1])
            if (dxHdz >= 0.246):
                dxHdz_xi_dec_smoothed_array.append(dxHdz)
                z_xi_dec_smoothed_array.append(antisym_func.dxHdz_To_z(dxHdz, M_max, zeta_z_func, T_vir, mu, z_dxHdz_history, dxHdz_history)[0])
    antisym_func.mkdir('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g'%(zeta, T_vir, R_mfp, BOX_LEN))
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g/dxHdz_xi_acc_smoothed_array.txt'%(zeta, T_vir, R_mfp, BOX_LEN), dxHdz_xi_acc_smoothed_array)
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g/dxHdz_xi_dec_smoothed_array.txt'%(zeta, T_vir, R_mfp, BOX_LEN), dxHdz_xi_dec_smoothed_array)
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g/z_xi_acc_smoothed_array.txt'%(zeta, T_vir, R_mfp, BOX_LEN), z_xi_acc_smoothed_array)
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g/z_xi_dec_smoothed_array.txt'%(zeta, T_vir, R_mfp, BOX_LEN), z_xi_dec_smoothed_array)
    
    #smooth the antisymmetric cross-correlation
    r12_smoothed_grid = np.logspace(-1, np.log10(300), 50)
    xi_acc_smoothed_func_array = []; xi_dec_smoothed_func_array = [] #[z](r12)
    for z in z_xi_acc_smoothed_array:
        xi_acc_smoothed_array = []
        for r12 in r12_smoothed_grid:
            xi_acc_smoothed_array.append(antisym_func.xi_A_HICO_smoothing(z, r12, xi_A_HICO_unsmoothed_func, BOX_LEN)[0])
        xi_acc_smoothed_func_array.append(interp1d(r12_smoothed_grid, xi_acc_smoothed_array, kind = 'cubic'))
    for z in z_xi_dec_smoothed_array:
        xi_dec_smoothed_array = []
        for r12 in r12_smoothed_grid:
            xi_dec_smoothed_array.append(antisym_func.xi_A_HICO_smoothing(z, r12, xi_A_HICO_unsmoothed_func, BOX_LEN)[0])
        xi_dec_smoothed_func_array.append(interp1d(r12_smoothed_grid, xi_dec_smoothed_array, kind = 'cubic'))

    #Fourier transformation to compute the smoothed antisymmetric power spectrum
    kh_array = np.logspace(np.log10(0.1),np.log10(0.6), 16)
    k_array = kh_array * antisym_func.hlittle
    Pk_A_acc = []; Pk_A_dec = [] #[z][kh]
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g/kh.txt'%(zeta, T_vir, R_mfp, BOX_LEN), kh_array, fmt='%9.9g')
    
    MP = Pool(NUM_CORE)
    for i in range(len(z_xi_acc_smoothed_array)):
        z = z_xi_acc_smoothed_array[i]
        MP.apply_async(Pk_A_function.Pk_A_multi_function, args = (k_array, xi_acc_smoothed_func_array[i], min(r12_smoothed_grid), max(r12_smoothed_grid), zeta, T_vir, R_mfp, z, BOX_LEN, ))
    MP.close()
    MP.join()
    
    MP = Pool(NUM_CORE)
    for i in range(len(z_xi_dec_smoothed_array)):
        z = z_xi_dec_smoothed_array[i]
        MP.apply_async(Pk_A_function.Pk_A_multi_function, args = (k_array, xi_dec_smoothed_func_array[i], min(r12_smoothed_grid), max(r12_smoothed_grid), zeta, T_vir, R_mfp, z, BOX_LEN, ))
    MP.close()
    MP.join()
    
    
    
    
