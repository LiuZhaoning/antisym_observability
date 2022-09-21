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
  print(z)
  tick1 = time.time()
  rhoHI_over_rho0 = antisym_func.HIrho_over_rho0(z, zeta_z_func, T_vir, mu, M_max)
  print(z, rhoHI_over_rho0, time.time() - tick1)
  results.put([z, rhoHI_over_rho0])
  

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
    #create an dir to restore the data
    DIR = '/scratch/liuzhaoning/antisym_observability/zeta%2.2g_Tvir%3.3g_Rmfp%2.2g_SMO%3.3g'%(zeta, T_vir, R_mfp, SMOOTHING_SCALE)
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
      z_HI = list(data['z_HI']); HI_ana = list(data['HI']); z_dxHdz = list(data['z_dxHdz']); dxHdz_ana = list(data['dxHdz'])
    else:
      z_HI = np.linspace(5.5, 12, 100); HI_ana = []
      for z in z_HI:
        HI_ana.append(1 - antisym_func.bar_Q(z, M_max, zeta_z_func, T_vir, mu, antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)))
      z_dxHdz, dxHdz_ana = antisym_func.dxH_dz_cal(z_HI, HI_ana)
      np.savez(DIR + '/history', z_HI = z_HI, HI = HI_ana, z_dxHdz = z_dxHdz, dxHdz = dxHdz_ana)

    #the redshift we are interested, mainly in the accelaration stage (dxHdz > 0.24)
    z_floor_xi_smoothed = z_dxHdz[dxHdz_ana.index(max(dxHdz_ana))]
    z_top_xi_smoothed = antisym_func.dxHdz_To_z(0.24, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana)[1]
    
    #the lowest redshift for xi_A_HICO_unsmoothing calculation considering the smoothing scale
    z_floor_xi_unsmoothed = antisym_func.cal_z1_z2(z_floor_xi_smoothed, SMOOTHING_SCALE, 0)[0] - 0.05 #room for conveniance
    z_top_xi_unsmoothed = antisym_func.cal_z1_z2(z_top_xi_smoothed, SMOOTHING_SCALE, 0)[1] + 0.05
    #calculate the redshift range considering both of smoothing and r12
    r12_limit = 150
    z_floor_HIrho = antisym_func.cal_z1_z2(z_floor_xi_unsmoothed, r12_limit, 0)[0] - 0.02
    z_top_HIrho = antisym_func.cal_z1_z2(z_top_xi_unsmoothed, r12_limit, 0)[1] + 0.02
    
    #compute anverage density of the neutral region
    if not os.path.exists(DIR + '/rhoHI_over_rho0.npz'):
      data = np.load(DIR + '/rhoHI_over_rho0.npz')
      z_array = data['z']; rhoHI_over_rho0_array = data['rhoHI']
    else:
      tick_start = time.time()
      z_array = np.linspace(z_floor_HIrho, z_top_HIrho, 5); rhoHI_over_rho0_array = np.zeros(100)
      MP = Pool(NUM_CORE); results = Manager().Queue()
      print('start')
      for z in z_array:
        MP.apply_async(HIrho_multiprocessing, args=(z, zeta_z_func, T_vir, R_mfp, mu, M_max, results,))
      MP.close(); MP.join()
      while not results.empty():
        z, rhoHI = results.get()
        print(z, rhoHI)
        rhoHI_over_rho0_array[z_array.index(z)] = rhoHI
      np.savez(DIR + '/rhoHI_over_rho0', z = z_array, rhoHI = rhoHI_over_rho0_array)
      print('the computation of average density of neutral region cost %3.3g mins'%((time.time() - tick_start) / 60))    
  



