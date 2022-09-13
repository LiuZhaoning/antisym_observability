#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from xi_A_HICO_unsmoothed_calculation_function import xi_A_HICO_z
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

#the parameters are choosed to match the history of the cyan line in Zhou et al. 2021

if __name__ == '__main__':
    zeta = float(sys.argv[1])
    T_vir = float(sys.argv[2])
    R_mfp = float(sys.argv[3])
    NUM_CORE = int(sys.argv[4])
    
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

    #read in the history of Zhou et al. (2021)
    z_HI = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/zhou/history/history_z.txt')
    HI_zhou = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/zhou/history/history_HI_cyan.txt')
    #calculate the dx_HI/dz
    [z_dxHdz, dxHdz_zhou] = antisym_func.dxH_dz_cal(z_HI, HI_zhou)

    #calculate the history
    HI_ana = []
    for z in z_HI:
        HI_ana.append(1 - antisym_func.bar_Q(z, M_max, zeta_z_func, T_vir, mu, antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)))
    [z_dxHdz, dxHdz_ana] = antisym_func.dxH_dz_cal(z_HI, HI_ana)
    print(HI_ana, dxHdz_ana)
    print('the maximum dxH/dz of the analytical model:', round(max(dxHdz_ana), 4))
    print('the maximum dxH/dz of 21cmFAST simulation:', round(max(dxHdz_zhou), 4))

    def dxHdz_To_z(dxHdz, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana):
        if (dxHdz >= max(dxHdz_ana) or dxHdz <= 0):
            print('The value should be an positive value smaller than the maximum dxHI/dz of the history')
            return 0
        z_turn = z_dxHdz[dxHdz_ana.index(max(dxHdz_ana))]
        
        def solution(z_floor, z_top):
            while (abs(z_floor - z_top) > 0.01):
                z_mid = (z_floor + z_top) / 2
                if (dxHdz > antisym_func.dxH_dz(z_mid, M_max, zeta_z_func, T_vir, mu)):
                    z_floor = z_mid
                else:
                    z_top = z_mid
            return (z_floor + z_top) / 2
        z1 = solution(6, z_turn)
        z2 = solution(10, z_turn)
        return [z1, z2]

    #the redshift of the lowest value of dxHI/dz in Zhou et al. 2021 (0.222 in acc and 0.246 in dec)
    z1 = dxHdz_To_z(0.222, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana)[0]
    z2 = dxHdz_To_z(0.246, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana)[1]
    #the lowest redshift for calculation considering smoothing in 384 Mpc box
    z_lower_limit = antisym_func.cal_z1_z2(z1, 384, 0)[0] - 0.05 #room for conveniance
    z_upper_limit = antisym_func.cal_z1_z2(z2, 384, 0)[1] + 0.05
    print('the redshift range for xi_A_HICO computation is (%4.4g, %4.4g)'%(z_lower_limit, z_upper_limit))
    #the minimum mass of bubbles in this redshift range
    M_min = zeta_z_func(z_upper_limit) * antisym_func.TtoM(z_upper_limit, T_vir, mu) * 0.99 #room for conveniance
    print('the minimum mass of the bubbles is %4.4g'%(M_min))

    #load in the BMF function
    z_grid_BMF = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g/z_grid_BMF.txt'%(zeta, T_vir, R_mfp))
    m_grid_BMF = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g/m_grid_BMF.txt'%(zeta, T_vir, R_mfp))
    dn_dm_mz = []
    for z in z_grid_BMF:
        dn_dm_mz.append(np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g/dec_BMF_grid_z%4.4g.txt'%(zeta, T_vir, R_mfp, z)))
    BMF_func = interp2d(m_grid_BMF, z_grid_BMF, dn_dm_mz, kind = 'cubic')

    #check the loaded data
    print('check the BMF_func')
    print(BMF_func(1e12, 7.3)[0], antisym_func.BMF(1e12, antisym_func.PARA_z(7.3, M_max, zeta_z_func, T_vir, mu)))

    #load in the HIrho_over_rho0 data
    z_array = np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/HIrho_over_rho0_zeta%3.3g_T%5.5g_R%2.2g/z_array.txt'%(zeta, T_vir, R_mfp))
    HIrho = []
    for i in range(len(z_array)):
        HIrho.append(np.loadtxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/HIrho_over_rho0_zeta%3.3g_T%5.5g_R%2.2g/HIrho_z%3.3g.txt'%(zeta, T_vir, R_mfp, z_array[i].flat[0])))

    #values at very low redshift is incredible
    for i in range(len(z_array)):
        if (HIrho[len(z_array) - 1 - i] < HIrho[len(z_array) - 2 - i]):
            for j in range(len(z_array) - 1 - i):
                HIrho[j] = HIrho[len(z_array) - 1 - i]

    HIrho_over_rho0_func = interp1d(z_array, HIrho, kind = 'cubic')
    #check the loaded data
    print('check the HIrho_over_rho0_func')
    print(HIrho_over_rho0_func(7.82), antisym_func.HIrho_over_rho0(7.82, zeta_z_func, T_vir, mu, M_max))
    
    #compute the antisymmetric cross correlation result at different redshift
    z_grid = np.linspace(z_lower_limit, z_upper_limit, 200)
    r12_limit = 150
    r12_grid = np.zeros(100)
    r12_grid[0:30] = np.linspace(0.1, 5, 30)
    r12_grid[30:100] = np.linspace(5, r12_limit, 71)[1:71]
    
    antisym_func.mkdir('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g'%(zeta, T_vir, R_mfp))
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g/z_grid'%(zeta, T_vir, R_mfp), z_grid, fmt = '%9.9g')
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g/r12_grid'%(zeta, T_vir, R_mfp), r12_grid, fmt = '%9.9g')
    
    tick_start = time.time()
    MP = Pool(NUM_CORE)
    for z in z_grid:
        MP.apply_async(xi_A_HICO_z, args=(z, zeta, zeta_z_func, HIrho_over_rho0_func, BMF_func, M_max, T_vir, R_mfp, mu, eta, r12_grid,))
    MP.close()
    MP.join()
    tick_end = time.time()
    print('xi_A_HICO_unsmoothed calculation cost %4.4g mins'%((tick_end - tick_start) / 60))
