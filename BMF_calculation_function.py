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
from multiprocessing import Pool
import antisym_func

#function for paralizing
def BMF_z(z, M_max, zeta, zeta_z_func, T_vir, R_mfp, mu, m_grid_BMF):
    tick_1 = time.time()
    dn_dm_array = []
    PARA = antisym_func.PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    for m in m_grid_BMF:
        dn_dm_array.append( antisym_func.BMF(m, PARA) )
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/BMF_zeta%3.3g_T%5.5g_R%2.2g/dec_BMF_grid_z%4.4g.txt'%(zeta, T_vir, R_mfp, z), dn_dm_array, fmt = '%.9g')
    tick_2 = time.time()
    print('BMF calculation at z=%5.5g cost time %3.3g min'%(z, (tick_2 - tick_1) / 60) )
    return 0
    
