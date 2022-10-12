import numpy as np
from scipy.interpolate import interp1d, interp2d
import time
import sys
import os
from multiprocessing import Pool, Manager
import antisym_func


def Pk_A_multiprocess(z, k_array, xi_A_smoothed_func, lower_limit, upper_limit, result):
    tick = time.time()
    Pk_A_array = []
    for k in k_array:
        Pk_A_array.append(antisym_func.Pk_A_cal(k, xi_A_smoothed_func, lower_limit, upper_limit))
    result.put([z, Pk_A_array])
    print('Pk_A compuatation at z=%3.3g cost %3.3g seconds'%(z, time.time() - tick))
    return 0

if __name__ == '__main__':
    tick_start = time.time()
    #parameters used in the computation
    zeta = float(sys.argv[1])
    T_vir = float(sys.argv[2])
    R_mfp = float(sys.argv[3])
    SMOOTHING_SCALE = float(sys.argv[4]) #Mpc
    SMOOTHING_Pk = float(sys.argv[5])
    NUM_CORE = int(sys.argv[6])
    M_max = antisym_func.RtoM(R_mfp)
    mu = 1.22 if T_vir < 9.99999e3 else 0.6
    
    #create an dir to restore the data
    antisym_func.mkdir('/scratch/liuzhaoning/antisym_observability/Pk_A_group_class')
    DIR_read = '/scratch/liuzhaoning/antisym_observability/xi_A_HICO/zeta%05.5g_Tvir%05.5g_Rmfp%05.5g_SMO%03.3g'%(zeta, T_vir, R_mfp, SMOOTHING_SCALE)
    DIR_load = '/scratch/liuzhaoning/antisym_observability/Pk_A_group_class/zeta%05.5g_Tvir%05.5g_Rmfp%05.5g_SMO%03.3g'%(zeta, T_vir, R_mfp, SMOOTHING_Pk)
    antisym_func.mkdir(DIR_load)
    
    #load in the normalized zeta
    if os.path.exists(DIR_read + '/normalized_zeta.npz'):
        data = np.load(DIR_read + '/normalized_zeta.npz')
        z_zeta_interp_array = data['z']; zeta_z_interp_array = data['zeta_z']
        zeta_z_func = interp1d(z_zeta_interp_array, zeta_z_interp_array, kind = 'cubic')
    else:
        print('file %s not found'%(DIR_read + '/nornalized_zeta.npz'))
        sys.exit()
    #load in the history
    if os.path.exists(DIR_read + '/history.npz'):
        data = np.load(DIR_read + '/history.npz')
        z_array_history = list(data['z_array_history']); HI_history = list(data['HI_history'])
        z_dxHdz_history = list(data['z_dxHdz_history']); dxHdz_history = list(data['dxHdz_history'])
        max_dxHdz = max(dxHdz_history)
    else:
        print('file %s not found'%(DIR_read + '/history.npz'))
        sys.exit()
    #load in the unsmoothed antisymmetric cross-correlation data
    if os.path.exists(DIR_read + '/xi_A_HICO_unsmoothed_map.npz'):
        data = np.load(DIR_read + '/xi_A_HICO_unsmoothed_map.npz')
        z_grid = list(data['z_grid']); r12_grid = list(data['r12_grid']); xi_A_HICO_map = list(data['xi_A_HICO_map'])       
        xi_A_HICO_unsmoothed_interp = interp2d(r12_grid, z_grid, xi_A_HICO_map, kind = 'cubic')
    else:
        print('file %s not found'%(DIR_read + '/xi_A_HICO_unsmoothed_map.npz'))
        sys.exit()
    
    #set the redshift of antisymmetric power spectrum we are computating
    dxHdz_temp_array = np.linspace(0.225, 0.445, 23)
    z_xi_acc_smoothed_array = []; z_xi_dec_smoothed_array = [];
    dxHdz_xi_acc_smoothed_array = []; dxHdz_xi_dec_smoothed_array = []
    for dxHdz in dxHdz_temp_array:
        if (dxHdz < max_dxHdz):
            dxHdz_xi_acc_smoothed_array.append(dxHdz); dxHdz_xi_dec_smoothed_array.append(dxHdz)
            z_xi_acc_smoothed_array.append(antisym_func.dxHdz_To_z(dxHdz, M_max, zeta_z_func, T_vir, mu, z_dxHdz_history, dxHdz_history)[1])
            z_xi_dec_smoothed_array.append(antisym_func.dxHdz_To_z(dxHdz, M_max, zeta_z_func, T_vir, mu, z_dxHdz_history, dxHdz_history)[0])
    
    #smooth the antisymmetric cross-correlation
    r12_limit = 150
    r12_smoothed_grid = np.zeros(100); r12_smoothed_grid[0:30] = np.linspace(0.1, 5, 30)
    r12_smoothed_grid[30:100] = np.linspace(5, r12_limit, 71)[1:71]
    xi_acc_smoothed_func_array = []; xi_dec_smoothed_func_array = [] #[z](r12)
    for z in z_xi_acc_smoothed_array:
        xi_acc_smoothed_array = []
        for r12 in r12_smoothed_grid:
            xi_acc_smoothed_array.append(antisym_func.xi_A_HICO_smoothing(z, r12, xi_A_HICO_unsmoothed_interp, SMOOTHING_Pk)[0])
        xi_acc_smoothed_func_array.append(interp1d(r12_smoothed_grid, xi_acc_smoothed_array, kind = 'cubic'))
    for z in z_xi_dec_smoothed_array:
        xi_dec_smoothed_array = []
        for r12 in r12_smoothed_grid:
            xi_dec_smoothed_array.append(antisym_func.xi_A_HICO_smoothing(z, r12, xi_A_HICO_unsmoothed_interp, SMOOTHING_Pk)[0])
        xi_dec_smoothed_func_array.append(interp1d(r12_smoothed_grid, xi_dec_smoothed_array, kind = 'cubic'))
  
    #Fourier transformation to compute the smoothed antisymmetric power spectrum
    kh_array = np.logspace(np.log10(0.1),np.log10(0.6), 16)
    k_array = kh_array * antisym_func.hlittle
    Pk_A_acc_map = [0] * len(z_xi_acc_smoothed_array); Pk_A_dec_map = [0] * len(z_xi_dec_smoothed_array) #[z][kh]
    
    MP = Pool(NUM_CORE); Pk_A_queue = Manager().Queue()
    for i in range(len(z_xi_acc_smoothed_array)):
        z = z_xi_acc_smoothed_array[i]
        MP.apply_async(Pk_A_multiprocess, args = (z, k_array, xi_acc_smoothed_func_array[i], min(r12_smoothed_grid), max(r12_smoothed_grid), Pk_A_queue ))
    MP.close(); MP.join()
    while not Pk_A_queue.empty():
        z, Pk_A_array = Pk_A_queue.get()
        Pk_A_acc_map[z_xi_acc_smoothed_array.index(z)] = Pk_A_array
    
    MP = Pool(NUM_CORE)
    for i in range(len(z_xi_dec_smoothed_array)):
        z = z_xi_dec_smoothed_array[i]
        MP.apply_async(Pk_A_multiprocess, args = (z, k_array, xi_dec_smoothed_func_array[i], min(r12_smoothed_grid), max(r12_smoothed_grid), Pk_A_queue ))
    MP.close(); MP.join()
    while not Pk_A_queue.empty():
        z, Pk_A_array = Pk_A_queue.get()
        Pk_A_dec_map[z_xi_dec_smoothed_array.index(z)] = Pk_A_array
    np.savez(DIR_load + '/Pk_A_acc_array', kh_array = kh_array, z_xi_acc_smoothed_array = z_xi_acc_smoothed_array, \
                dxHdz_xi_acc_smoothed_array = dxHdz_xi_acc_smoothed_array, Pk_A_acc_map = Pk_A_acc_map)
    np.savez(DIR_load + '/Pk_A_dec_array', z_xi_dec_smoothed_array = z_xi_dec_smoothed_array, \
                dxHdz_xi_dec_smoothed_array = dxHdz_xi_dec_smoothed_array, Pk_A_dec_map = Pk_A_dec_map)
                
    print('Pk_A_HICO calculation cost %4.4g mins'%((time.time() - tick_start) / 60))
