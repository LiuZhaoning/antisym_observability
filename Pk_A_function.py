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

def Pk_A_multi_function(k_array, xi_A_smoothed_func, lower_limit, upper_limit, zeta, T_vir, R_mfp, z, BOX_LEN):
    tick1 = time.time()
    Pk_A_array = []
    for k in k_array:
        Pk_A_array.append(antisym_func.Pk_A_cal(k, xi_A_smoothed_func, lower_limit, upper_limit))
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/Pk_A_HICO_zeta%3.3g_T%5.5g_R%2.2g_LEN%3.3g/Pk_A_z%3.3g.txt'%(zeta, T_vir, R_mfp, BOX_LEN, z), Pk_A_array, fmt = '%9.9g')
    tick2 = time.time()
    print('Pk_A compuatation at z=%3.3g cost %3.3g seconds'%(z, tick2 - tick1))
    return 0
    
    
    
