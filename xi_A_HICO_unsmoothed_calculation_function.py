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

def xi_A_HICO_z(z, zeta, zeta_z_func, HIrho_over_rho0_func, BMF_func, M_max, T_vir, R_mfp, mu, eta, r12_grid):
    tick1 = time.time()
    xi_A_HICO_array = []
    for r12 in r12_grid:
        xi_A_HICO_array.append(antisym_func.xi_A_HICO(z, r12, zeta_z_func, HIrho_over_rho0_func, BMF_func, M_max, T_vir, mu, eta))
    np.savetxt('/home/liuzhaoning/antisym/antisym_parameters/cal_mid_data/xi_A_HICO_unsmoothed_zeta%3.3g_T%5.5g_R%2.2g/xi_A_HICO_unsmoothed_z%4.4g.txt'%(zeta, T_vir, R_mfp, z), xi_A_HICO_array, fmt='%9.9g')
    tick2 = time.time()
    print('xi_A_HICO computation at z = %3.3g cost %3.3g mins'%(z, (tick2-tick1) / 60))
    return 0
