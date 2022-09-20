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
import antisym_func
from BMF_calculation_function import BMF_z

if __name__ == '__main__':
    #put in parameters of the reionization model
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
    NUM = 50; z_zeta_interp_array = np.linspace(5.5, 13.5, NUM); zeta_z_interp_array = []
    for i in range(NUM):
        zeta_z_interp_array.append(antisym_func.zeta_z(z_zeta_interp_array[i], zeta, T_vir, mu))
    #zeta(z) interpolation
    zeta_z_func = interp1d(z_zeta_interp_array, zeta_z_interp_array, kind = 'cubic')

    #calculate the history
    HI_ana = []
    for z in z_HI:
        HI_ana.append(1 - antisym_func.bar_Q(z, M_max, zeta_z_func, T_vir, mu, antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)))
    z_dxHdz, dxHdz_ana = antisym_func.dxH_dz_cal(z_HI, HI_ana)

    #the redshift we are interested, mainly in the accelaration stage (dxHdz > 0.24)
    z1 = z_dxHdz[dxHdz_ana.index(max(dxHdz_ana))] 
    z2 = dxHdz_To_z(0.24, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana)[1]
    #the lowest redshift for xi_A_HICO_unsmoothing calculation considering the smoothing scale
    z_lower_limit = antisym_func.cal_z1_z2(z1, SMOOTHING_SCALE, 0)[0] - 0.05 #room for conveniance
    z_upper_limit = antisym_func.cal_z1_z2(z2, SMOOTHING_SCALE, 0)[1] + 0.05
    #the minimum mass of bubbles in this redshift range
    M_min = zeta_z_func(z_upper_limit) * antisym_func.TtoM(z_upper_limit, T_vir, mu) * 0.99 #room for conveniance
    #calculate the redshift range considering both of smoothing and r12
    r12_limit = 150
    z_floor_HIrho = antisym_func.cal_z1_z2(z_lower_limit, r12_limit, 0)[0] - 0.02
    z_top_HIrho = antisym_func.cal_z1_z2(z_upper_limit, r12_limit, 0)[1] + 0.02
    print('the range of redshift of HIrho computation is (%3.3g, %3.3g)'%(z_floor_BMF, z_top_BMF))
    
    

    #calculate the dn/dm(z,m)
    MP = Pool(NUM_CORE)
    tick_start = time.time()
    for z in z_grid_BMF:
        MP.apply_async(BMF_z, args = (z, M_max, zeta, zeta_z_func, T_vir, R_mfp, mu, m_grid_BMF,))
    MP.close()
    MP.join()
    tick_end = time.time()
    print('BMF calculation cost time %3.3g min'%((tick_end - tick_start) / 60) )




