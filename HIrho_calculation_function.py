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

def HIrho_multiprocessing(z, zeta, zeta_z_func, T_vir, R_mfp, mu, M_max):
    tick1 = time.time()
    HIrho=[]
    HIrho.append(antisym_func.HIrho_over_rho0(z, zeta_z_func, T_vir, mu, M_max))
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/HIrho_over_rho0_zeta%3.3g_T%5.5g_R%2.2g/HIrho_z%3.3g.txt'%(zeta, T_vir, R_mfp, z), HIrho, fmt='%9.9g')
    tick2 = time.time()
    print('(1 + bar{delta}) = %3.3g at z=%3.3g, cost time %3.3g seconds'%(HIrho[0], z, tick2 - tick1))
    return 0

