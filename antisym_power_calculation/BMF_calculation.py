import numpy as np
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
from multiprocessing import Pool
import ../antisym_func
from BMF_calculation_function import BMF_z
#test
if __name__ == '__main__':
    #put in the number of cores for multiprocessing
    zeta = float(sys.argv[1])
    T_vir = float(sys.argv[2])
    R_mfp = float(sys.argv[3])
    NUM_CORE = int(sys.argv[4])

    #parameters used in the computation
    eta = 0.5 #the limit where we overlook the correlation outside bubbles
    M_max = antisym_func.RtoM(R_mfp)
    mu = 1.22 if T_vir < 9.99999e3 else 0.6
    SMOOTHING_SCALE = 300 #Mpc

    #calculate some of the zeta_z points to do interpolation
    NUM = 50
    z_zeta_interp_array = np.linspace(5.5, 13.5, NUM)
    zeta_z_interp_array = []
    for i in range(NUM):
        zeta_z_interp_array.append(antisym_func.zeta_z(z_zeta_interp_array[i], zeta, T_vir, mu))
    #zeta(z) interpolation
    zeta_z_func = interp1d(z_zeta_interp_array, zeta_z_interp_array, kind = 'cubic')

    #calculate the history
    HI_ana = []
    for z in z_HI:
        HI_ana.append(1 - antisym_func.bar_Q(z, M_max, zeta_z_func, T_vir, mu, antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)))
    [z_dxHdz, dxHdz_ana] = antisym_func.dxH_dz_cal(z_HI, HI_ana)

    #the redshift we are interested, mainly in the accelaration stage (dxHdz > 0.25)
    z1 = z_dxHdz[dxHdz_ana.index(max(dxHdz_ana))] 
    z2 = dxHdz_To_z(0.246, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana)[1]
    #the lowest redshift for xi_A_HICO_unsmoothing calculation considering the smoothing scale
    z_lower_limit = antisym_func.cal_z1_z2(z1, SMOOTHING_SCALE, 0)[0] - 0.05 #room for conveniance
    z_upper_limit = antisym_func.cal_z1_z2(z2, SMOOTHING_SCALE, 0)[1] + 0.05
    print('the redshift range for xi_A_HICO computation is (%4.4g, %4.4g)'%(z_lower_limit, z_upper_limit))
    #the minimum mass of bubbles in this redshift range
    M_min = zeta_z_func(z_upper_limit) * antisym_func.TtoM(z_upper_limit, T_vir, mu) * 0.99 #room for conveniance
    print('the minimum mass of the bubbles is %4.4g'%(M_min))
    
    #calculate the redshift range considering the r12, other than the smoothing
    r12_limit = 150
    z_floor_BMF = antisym_func.cal_z1_z2(z_lower_limit, r12_limit, 0)[0] - 0.02
    z_top_BMF = antisym_func.cal_z1_z2(z_upper_limit, r12_limit, 0)[1] + 0.02
    print('the range of redshift of BMF computation is (%3.3g, %3.3g)'%(z_floor_BMF, z_top_BMF))
    
    #calculate the dn/dm(z,m), the 2D array for interpolation
    #set up the 2D mesh index
    z_grid_BMF = np.linspace(z_floor_BMF, z_top_BMF, 500)
    m_grid_BMF = np.logspace(np.log10(M_min), np.log10(0.5 * M_max), 400) #m>6e8
    m_grid_BMF = np.resize(m_grid_BMF,500)
    m_grid_BMF[400:500] = np.linspace(0.5 * M_max, 0.999 * M_max, 101)[1:101]

    #save the grid index
    antisym_func.mkdir('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g'%(zeta, T_vir, R_mfp))
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g/z_grid_BMF.txt'%(zeta, T_vir, R_mfp), z_grid_BMF, fmt = '%.9g')
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g/m_grid_BMF.txt'%(zeta, T_vir, R_mfp), m_grid_BMF, fmt = '%.9g')

    #calculate the dn/dm(z,m)
    MP = Pool(NUM_CORE)
    tick_start = time.time()
    for z in z_grid_BMF:
        MP.apply_async(BMF_z, args = (z, M_max, zeta, zeta_z_func, T_vir, R_mfp, mu, m_grid_BMF,))
    MP.close()
    MP.join()
    tick_end = time.time()
    print('BMF calculation cost time %3.3g min'%((tick_end - tick_start) / 60) )




