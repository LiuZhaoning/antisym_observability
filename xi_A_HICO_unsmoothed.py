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
from multiprocessing import Pool, Manager
import antisym_func

def HIrho_multiprocessing(z, zeta_z_func, T_vir, R_mfp, mu, M_max, results):
  tick1 = time.time()
  rhoHI_over_rho0 = antisym_func.HIrho_over_rho0(z, zeta_z_func, T_vir, mu, M_max)
  results.put([z, rhoHI_over_rho0])
  print(z, 'finish HIrho cost', (time.time() - tick1) / 60, 'mins')

def xi_A_HICO_z(z, zeta_z_func, HIrho_over_rho0_func, M_max, T_vir, R_mfp, mu, eta, r12_grid, m_array_BMF, results):
    tick1 = time.time()
    BMF_array = []
    PARA = antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    for m in m_array_BMF:
        BMF_array.append(antisym_func.BMF(m, PARA))
    BMF_func_z = interp1d(m_array_BMF, BMF_array, kind = 'cubic')
    xi_A_HICO_array = []
    for r12 in r12_grid:
        xi_A_HICO_array.append(antisym_func.xi_A_HICO(z, r12, zeta_z_func, HIrho_over_rho0_func, BMF_func_z, M_max, T_vir, mu, eta))
    results.put([z, BMF_array, xi_A_HICO_array])
    print(z, 'finish xi_A_HICO cost', (time.time() - tick1) / 60, 'mins')

if __name__ == '__main__':
    #put in parameters of the reionization model
    zeta = float(sys.argv[1])
    T_vir = float(sys.argv[2])
    R_mfp = float(sys.argv[3])
    SMOOTHING_SCALE = float(sys.argv[4]) #Mpc
    NUM_CORE = int(sys.argv[5])

    #parameters used in the computation
    eta = 0.5 #the limit where we overlook the correlation outside bubbles
    M_max = antisym_func.RtoM(R_mfp)
    mu = 1.22 if T_vir < 9.99999e3 else 0.6
    #create an dir to restore the data
    DIR = '/scratch/liuzhaoning/antisym_observability/xi_A_HICO/zeta%05.5g_Tvir%05.5g_Rmfp%05.5g_SMO%03.3g'%(zeta, T_vir, R_mfp, SMOOTHING_SCALE)
    #check whether the condition is already computated
    if os.path.exists(DIR + '/xi_A_HICO_unsmoothed_map.npz'): 
      print('the condition zeta%5.5g_Tvir%5.5g_Rmfp%5.5g_SMO%3.3g is already computated'(zeta, T_vir, R_mfp, SMOOTHING_SCALE)) 
    antisym_func.mkdir(DIR)

    #calculate normalized zeta_z points to do interpolation
    if os.path.exists(DIR + '/normalized_zeta.npz'): 
      data = np.load(DIR + '/normalized_zeta.npz')
      z_zeta_interp_array = data['z']; zeta_z_interp_array = data['zeta_z']
    else:
      tick_start = time.time()
      NUM = 100; z_zeta_interp_array = np.linspace(5.5, 12, NUM); zeta_z_interp_array = []
      for i in range(NUM):
        zeta_z_interp_array.append(antisym_func.zeta_z(z_zeta_interp_array[i], zeta, T_vir, mu))
      np.savez(DIR + '/normalized_zeta', z = z_zeta_interp_array, zeta_z = zeta_z_interp_array)
      print('the computation of normalized zeta cost %3.3g mins'%((time.time() - tick_start) / 60))
    #zeta(z) interpolation
    zeta_z_func = interp1d(z_zeta_interp_array, zeta_z_interp_array, kind = 'cubic')

    #calculate the history of HI fraction and dx_HI/dz
    if os.path.exists(DIR + '/history.npz'):
        data = np.load(DIR + '/history.npz')
        z_array_history = list(data['z_array_history']); HI_history = list(data['HI_history'])
        z_dxHdz_history = list(data['z_dxHdz_history']); dxHdz_history = list(data['dxHdz_history'])
    else:
        z_array_history = np.linspace(5.5, 12, 100); HI_history = []
        for z in z_array_history:
            HI_history.append(1 - antisym_func.bar_Q(z, M_max, zeta_z_func, T_vir, mu, antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)))
        z_dxHdz_history, dxHdz_history = antisym_func.dxH_dz_cal(z_array_history, HI_history)
        np.savez(DIR + '/history', z_array_history = z_array_history, HI_history = HI_history, \
                    z_dxHdz_history = z_dxHdz_history, dxHdz_history = dxHdz_history)

    #the redshift we are interested, mainly in the accelaration stage (dxHdz > 0.24)
    z_floor_xi_smoothed = antisym_func.dxHdz_To_z(0.222, M_max, zeta_z_func, T_vir, mu, z_dxHdz_history, dxHdz_history)[0]
    z_top_xi_smoothed = antisym_func.dxHdz_To_z(0.222, M_max, zeta_z_func, T_vir, mu, z_dxHdz_history, dxHdz_history)[1]
    
    #the lowest redshift for xi_A_HICO_unsmoothing calculation considering the smoothing scale
    z_floor_xi_unsmoothed = antisym_func.cal_z1_z2(z_floor_xi_smoothed, SMOOTHING_SCALE, 0)[0] - 0.05 #room for conveniance
    z_top_xi_unsmoothed = antisym_func.cal_z1_z2(z_top_xi_smoothed, SMOOTHING_SCALE, 0)[1] + 0.05
    #the minimum mass of bubbles in this redshift range
    M_min = zeta_z_func(z_top_xi_unsmoothed) * antisym_func.TtoM(z_top_xi_unsmoothed, T_vir, mu) * 0.99 #room for conveniance
    #calculate the redshift range considering both of smoothing and r12
    r12_limit = 150
    z_floor_HIrho = antisym_func.cal_z1_z2(z_floor_xi_unsmoothed, r12_limit, 0)[0] - 0.02
    z_top_HIrho = antisym_func.cal_z1_z2(z_top_xi_unsmoothed, r12_limit, 0)[1] + 0.02
    
    #compute anverage density of the neutral region
    if os.path.exists(DIR + '/rhoHI_over_rho0.npz'):
        data = np.load(DIR + '/rhoHI_over_rho0.npz')
        z_array_HIrho = list(data['z_array_HIrho']); rhoHI_over_rho0_array = list(data['rhoHI_over_rho0_array'])
    else:
        NUM = 100; tick_start = time.time()
        z_array_HIrho = list(np.linspace(z_floor_HIrho, z_top_HIrho, NUM)); rhoHI_over_rho0_array = list(np.zeros(NUM))
        MP = Pool(NUM_CORE); rhoHI_queue = Manager().Queue()
        print('start')
        for z in z_array_HIrho:
            MP.apply_async(HIrho_multiprocessing, args=(z, zeta_z_func, T_vir, R_mfp, mu, M_max, rhoHI_queue,))
        MP.close(); MP.join()
        while not rhoHI_queue.empty():
            z, rhoHI = rhoHI_queue.get()
            rhoHI_over_rho0_array[z_array_HIrho.index(z)] = rhoHI
        np.savez(DIR + '/rhoHI_over_rho0', z_array_HIrho = z_array_HIrho, rhoHI_over_rho0_array = rhoHI_over_rho0_array)
        print('the computation of average density of neutral region cost %3.3g mins'%((time.time() - tick_start) / 60))
    #values at very low redshift is incredible, so I set them to be the minimum HIrho_over_rho0
    for i in range(len(z_array_HIrho)):
        if (rhoHI_over_rho0_array[len(z_array_HIrho) - 1 - i] < rhoHI_over_rho0_array[len(z_array_HIrho) - 2 - i]):
            for j in range(len(z_array_HIrho) - 1 - i):
                rhoHI_over_rho0_array[j] = rhoHI_over_rho0_array[len(z_array_HIrho) - 1 - i]
    HIrho_over_rho0_interp = interp1d(z_array_HIrho, rhoHI_over_rho0_array, kind = 'cubic')

    #use the m_array to calculate the BMF(m, const z) interpolation
    m_array_BMF = np.logspace(np.log10(M_min), np.log10(0.5 * M_max), 200) 
    m_array_BMF = np.resize(m_array_BMF,250)
    m_array_BMF[200:250] = np.linspace(0.5 * M_max, 0.999 * M_max, 51)[1:51]
    #set up the grid for z and r12
    z_grid = list(np.linspace(z_floor_xi_unsmoothed, z_top_xi_unsmoothed, 200))
    r12_grid = np.zeros(100); r12_grid[0:30] = np.linspace(0.1, 5, 30); r12_grid[30:100] = np.linspace(5, r12_limit, 71)[1:71]
    BMF_map = [0] * len(z_grid); xi_A_HICO_map = [0] * len(z_grid)

    #computate the xi_A_HICO grid
    tick_start = time.time()
    MP = Pool(NUM_CORE); xi_queue = Manager().Queue()
    for z in z_grid:
        MP.apply_async(xi_A_HICO_z, args=(z, zeta_z_func, HIrho_over_rho0_interp, M_max, \
                                          T_vir, R_mfp, mu, eta, r12_grid, m_array_BMF, xi_queue,))
    MP.close(); MP.join()
    while not xi_queue.empty():
        z, BMF_array, xi_A_HICO_array = xi_queue.get()
        index = z_grid.index(z)
        BMF_map[index] = BMF_array; xi_A_HICO_map[index] = xi_A_HICO_array
    np.savez(DIR + '/BMF_map', z_grid = z_grid, m_grid = m_array_BMF, BMF_map = BMF_map)
    np.savez(DIR + '/xi_A_HICO_unsmoothed_map', z_grid = z_grid, r12_grid = r12_grid, xi_A_HICO_map = xi_A_HICO_map)
    print('xi_A_HICO_unsmoothed calculation cost %4.4g mins'%((time.time() - tick_start) / 60))
