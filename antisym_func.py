import numpy as np
import math
from scipy import integrate
from scipy import special
from scipy.interpolate import interp1d, interp2d
from scipy.misc import derivative
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import time
import mpmath
import os

#change

speed_C = 2.99792e5 # km/s
hlittle = (float)(0.678)
SIGMA8 = 0.82
P_CUTOFF = 0
FILTER = 0
POWER_INDEX = 0.968
SHETH_a = (float)(0.73)
SHETH_p = (float)(0.173)
SHETH_A = (float)(0.353)
E = (float)(2.71828182846)
OMn = (float)(0.0)
OMm = (float)(0.306)
OMl = (float)(1.0-OMm)
OMr = (float)(8.6e-5)
OMb = (float) (0.0492)
OMk = 0.0
TINY = (float)(1e-30)
T_cmb = 2.728
N_nu = 1.0
g_x = 1.5
M_WDM = 2.0

BODE_e = 0.361
BODE_n = 5.0
BODE_v = 1.2

Unit_SolarMass = 1.9891e30 #kg
Unit_Mpc = 3.08567758e22 #m
const_k = 1.380649e-23 # kg*m^2*s^-2
G = 6.6743015e-11 # kg^-1*m^3*s^-2
m_H = 1.6735575e-27 #kg
#solar_mass = 1.9891e33 #g
#Mpc = 3.08567758e25 #cm
rho_crit_0 = 2.7754e11 * pow(hlittle,2) #solarmass per Mpc^-3
delta_c = 1.68
rho_bar = rho_crit_0 * OMm # solarmass per Mpc^-3 #for matter

#calculate the P^I_HICO with multiprocess
def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
        print ("---  new folder created  ---")
    else:
        print ("---  The folder has existed  ---")
    return 0

#calculate the power spectrum at z=0 in the same way with 21cmFAST (POWERSPECTRM = 0)
R_CUTOFF = 0.201*pow((OMm-OMb)*hlittle*hlittle/0.15, 0.15)*pow(g_x/1.5, -0.29)*pow(M_WDM, -1.15);
omhh = OMm*hlittle*hlittle;
theta_cmb = T_cmb / 2.7
f_nu = OMn/OMm;
f_baryon = OMb/OMm;

if (f_nu < TINY):
    f_nu = 1e-10
    
if (f_baryon < TINY):
    f_baryon = 1e-10
    
#TFset_Parameters in 21cmFAST
z_equality = 25000*omhh*pow(theta_cmb, -4) - 1.0;
k_equality = 0.0746*omhh/(theta_cmb*theta_cmb);

z_drag = 0.313*pow(omhh,-0.419) * (1 + 0.607*pow(omhh, 0.674));
z_drag = 1 + z_drag*pow(OMb*hlittle*hlittle, 0.238*pow(omhh, 0.223));
z_drag *= 1291 * pow(omhh, 0.251) / (1 + 0.659*pow(omhh, 0.828));

y_d = (1 + z_equality) / (1.0 + z_drag);

R_drag = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_drag);
R_equality = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_equality);

sound_horizon = 2.0/3.0/k_equality * math.sqrt(6.0/R_equality) * math.log( (math.sqrt(1+R_drag) + math.sqrt(R_drag+R_equality)) / (1.0 + math.sqrt(R_equality)) );
p_c = -(5 - math.sqrt(1 + 24*(1 - f_nu-f_baryon)))/4.0;
p_cb = -(5 - math.sqrt(1 + 24*(1 - f_nu)))/4.0;
f_c = 1 - f_nu - f_baryon;
f_cb = 1 - f_nu;
f_nub = f_nu+f_baryon;

alpha_nu = (f_c/f_cb) * (2*(p_c+p_cb)+5)/(4*p_cb+5.0);
alpha_nu *= 1 - 0.553*f_nub+0.126*pow(f_nub,3);
alpha_nu /= 1-0.193 * math.sqrt(f_nu)+0.169*f_nu;
alpha_nu *= pow(1+y_d, p_c-p_cb);
alpha_nu *= 1+ (p_cb-p_c)/2.0 * (1.0+1.0/(4.0*p_c+3.0)/(4.0*p_cb+7.0))/(1.0+y_d);
beta_c = 1.0/(1.0-0.949*f_nub);

sigma_norm = -1;

R = 8.0/hlittle;
kstart = 1.0e-99/R;
kend = 350.0/R;
lower_limit = kstart;#log(kstart);
upper_limit = kend;#log(kend);

def TFmdm(k):
    q = k * pow(theta_cmb,2)/omhh
    gamma_eff = math.sqrt(alpha_nu) + (1.0 - math.sqrt(alpha_nu)) / (1.0+pow(0.43 * k * sound_horizon, 4))    
    q_eff = q / gamma_eff
    TF_m= math.log(E+1.84*beta_c * math.sqrt(alpha_nu)*q_eff)
    TF_m /= TF_m + pow(q_eff,2) * (14.4 + 325.0/(1.0+60.5*pow(q_eff,1.11)))   
    q_nu = 3.92*q/ math.sqrt(f_nu/N_nu)
    TF_m *= 1.0 + (1.2*pow(f_nu,0.64)*pow(N_nu,0.3+0.6*f_nu)) / (pow(q_nu,-1.6)+pow(q_nu,0.8))
    return TF_m;

def dsigma_dk(k,R):
    T = TFmdm(k);
    # check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF):
        T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v)
    p = pow(k, POWER_INDEX) * T * T;
    kR = k*R;
    if ( (FILTER == 0) or (sigma_norm < 0) ):# top hat
        if ( (kR) < 1.0e-4 ):
            w = 1.0 # w converges to 1 as (kR) -> 0
        else:
            w = 3.0 * (np.sin(kR)/pow(kR, 3) - np.cos(kR)/pow(kR, 2))
    if (FILTER == 1):
        w = pow(E, -kR*kR/2.0)
    return k*k*p*w*w

result = integrate.quad(dsigma_dk, lower_limit, upper_limit,(R))
sigma_norm = SIGMA8 / (result[0] ** 0.5)

def power_in_k(k):#the power spectrum of matter density at z=0
    T = TFmdm(k)
    if (P_CUTOFF):
        T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v)
    p = pow(k, POWER_INDEX) * T * T;
    return p * (2 * math.pi * math.pi) * (sigma_norm * sigma_norm)
vec_power = np.vectorize(power_in_k)

#extend the power spectrum to the small scale
def largek_extention(k, A, alpha): #xi_extend(r) = A * r ^ alpha
    return A * k ** alpha
vec_larger_extention = np.vectorize(largek_extention)

#link the sigma to R and M
def MtoR(M):
    return pow(3 * M / (4 * math.pi * rho_bar), 1.0 / 3.0)
vec_MtoR = np.vectorize(MtoR)

#the mass of bubbles with radius R and average matter density
def RtoM(R):
    M = 4.0 / 3.0 * math.pi * pow(R,3) * rho_bar
    return M
M_max = RtoM(50) #default

#the variance square of the overdensity of bubbles with mass M 
def S(M):
    R = (M / (6 * np.pi ** 2 * rho_bar)) ** (1.0 / 3) #this R is only for estimation of k-space cutoff
    def integrand(k):
        return power_in_k(k) * k * k;
    sig2 = 1 / (2 * math.pi * math.pi) * integrate.quad(integrand, 0, 1 / R, epsrel=1e-5)[0] #I set a upper limit since the calculation of inf cost too much time
    return sig2

#the variance of the overdensity of bubbles with mass M 
def sigma(M):
    return S(M) ** 0.5
    
def omega_mz(z):
    return OMm*pow(1+z,3) / (OMm*pow(1+z,3) + OMl + OMr*pow(1+z,4) + OMk*pow(1+z, 2))

def Deltac_nonlinear(z):
    d = omega_mz(z) -1.0
    return 18 * np.pi * np.pi - 39 * d * d

#the virial temperature to the minimum mass of ionizing halos
def TtoM(z, T, mu): #in sun_mass
    return 7030.97 / hlittle * ( omega_mz(z) / (OMm*Deltac_nonlinear(z))) ** 0.5 * pow( T/(mu * (1+z)),1.5)

#calculate H(z)
def H_z(z):
    H_0 = 100 * hlittle # km/s/Mpc
    H_z = H_0 * pow(OMr * pow(1 + z,4) + OMm * pow(1 + z,3) + (1 - OMr - OMm) + OMk * pow(1 + z,2),0.5)
    return H_z;

#growth function in 21cmFAST
def dicke(z):
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) )
    dick_z = 2.5*omegaM_z / ( 1.0/70.0 + omegaM_z*(209-omegaM_z)/140.0 + pow(omegaM_z, 4.0/7.0) )
    dick_0 = 2.5*OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) )
    return dick_z / (dick_0 * (1.0+z))

#locate the different redshift
def LOS_distance(z): #the comoving distance from now to the redshift z
    def dchi(z): # dchi = c/H(z) dz
        return speed_C / H_z(z)
    chi = integrate.quad(dchi, 0, z, epsrel = 1e-5)
    return chi[0]

#at redshift z, the function will give the 2 redshift with comoving distance r/2
def cal_z1_z2(z, r12, theta):
    def functioin_z1(z1):
        return (LOS_distance(z) - LOS_distance(z1) - r12 * np.cos(theta) / 2)
    def functioin_z2(z2):
        return (LOS_distance(z) - LOS_distance(z2) + r12 * np.cos(theta) / 2)
    z1 = fsolve(functioin_z1, [z])[0]
    z2 = fsolve(functioin_z2, [z])[0]
    return [z1, z2]

#f_coll comes from ST distribution
def sigma_z0(M):
    R = pow(3 * M / (4 * math.pi * rho_bar),1.0/3.0)
    kstart = 1.0e-99/R
    kend = 350.0/R
    lower_limit = kstart
    upper_limit = kend
    result = integrate.quad(dsigma_dk, lower_limit, upper_limit, epsrel = 1e-5, args=(R))
    return sigma_norm * (result[0] ** 0.5)

def dsigmasq_dm(k, d2fact, R):
    T = TFmdm(k)
    if (P_CUTOFF):
        T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v)
    p = pow(k, POWER_INDEX) * T * T
    kR = k * R
    if (FILTER == 0):# top hat
        if ( (kR) < 1.0e-4 ):
            w = 1.0 # w converges to 1 as (kR) -> 0
        else:
            w = 3.0 * (np.sin(kR)/pow(kR, 3) - np.cos(kR)/pow(kR, 2))
        if ( (kR) < 1e-10 ):
            dwdr = 0
        else:
            dwdr = 9 * np.cos(kR) * k/pow(kR,3) + 3*np.sin(kR)*(1 - 3/(kR*kR))/(kR*R)
        drdm = 1.0 / (4.0 * np.pi * rho_bar * R * R)
    return k*k*p*2*w*dwdr*drdm * d2fact
    
def dsigmasqdm_z0(M):
    R = pow(3 * M / (4 * math.pi * rho_bar),1.0/3.0)
    kstart = 1.0e-99/R
    kend = 350.0/R
    lower_limit = kstart
    upper_limit = kend
    d2fact = M*10000/sigma_z0(M)
    result = integrate.quad(dsigmasq_dm, lower_limit, upper_limit,epsrel = 1e-3, args = (d2fact,R), limit=100)
    return sigma_norm * sigma_norm * result[0] / d2fact

def dNdM_st(z, M):
    dicke_growth = dicke(z);
    sigma = sigma_z0(M) * dicke_growth;
    dsigmadm = dsigmasqdm_z0(M) * dicke_growth*dicke_growth/(2.0*sigma);
    nuhat = (SHETH_a ** 0.5) * delta_c / sigma;
    return (-rho_bar/M) * (dsigmadm/sigma) * ((2./np.pi) ** 0.5) * SHETH_A * (1+ pow(nuhat, -2*SHETH_p)) * nuhat * pow(E, -nuhat*nuhat/2.0)

def dFdlnM_st(lnM,z):
    M = np.exp(lnM)
    result = dNdM_st(z, M) * M * M
    return result
    
def FgtrM_st(z,M):
    lower_limit = np.log(M)
    upper_limit = np.log(max(1e16, M * 100))
    result = integrate.quad(dFdlnM_st, lower_limit, upper_limit, epsrel = 1e-3,args = (z))
    return result[0] / (rho_bar)

#the collapsed fraction in ST method
def f_coll_st(z, T, mu):
    M_min = TtoM(z, T, mu)
    f_coll = FgtrM_st(z, M_min)
    return f_coll
vec_f_coll_st = np.vectorize(f_coll_st)

def X_H_fcoll_st(z,zeta,T):
    return (1 - zeta * f_coll_st(z,T));
vec_XH_fcoll_st = np.vectorize(X_H_fcoll_st)

def dXH_dz_fcoll_st(z,zeta,T):
    def XH(z): #for derivative only
        return X_H_fcoll_st(z,zeta,T)
    deriv = derivative(XH,z,dx=1e-2)
    return deriv
vec_dXH_dz_fcoll_st = np.vectorize(dXH_dz_fcoll_st)

#f_coll comes from PS distribution
def f_coll(z, T, mu):
    delta_c_z = delta_c / dicke(z)
    M_min = TtoM(z, T, mu)  
    f_coll = special.erfc(delta_c_z / (pow(2,0.5) * sigma(M_min)))
    return f_coll
vec_f_coll = np.vectorize(f_coll)

#since the value in 21cmFAST was normalized to f_coll_st, we consider zeta(z) = zeta * f_coll_st / f_coll
def zeta_z(z, zeta, T_vir, mu):
    return zeta * f_coll_st(z,T_vir, mu) / f_coll(z, T_vir, mu)

#Set the barrier
def delta_x(z, S0, zeta_z_func, T_vir, mu):
    delta_c_z = delta_c / dicke(z)
    M_min = TtoM(z, T_vir, mu)
    delta_x = delta_c_z - pow(2.0,0.5) * special.erfinv( 1.0 - 1.0 / zeta_z_func(z) ) * pow(S(M_min) - S0,0.5);
    return delta_x;

#set the fitted barrier
def B_0(z, zeta_z_func, T_vir, mu):
    delta_c_z = delta_c / dicke(z)
    M_min = TtoM(z, T_vir, mu)
    B_0 = delta_c_z - pow(2.0,0.5) * special.erfinv( 1.0 - 1.0 / zeta_z_func(z) ) * pow(S(M_min),0.5)
    return B_0;

def B_1(z, zeta_z_func, T_vir, mu):
    M_min = TtoM(z, T_vir, mu)
    slope = special.erfinv( 1.0 - 1.0 / zeta_z_func(z) ) / pow(2*S(M_min),0.5)
    return slope;

#set the curves end at S(zeta*M_min)
def Cutoff(z, zeta_z_func, T_vir, mu):
    M_min = TtoM(z, T_vir, mu)
    cut = S(zeta_z_func(z) * M_min)
    return cut;
vec_delta_x = np.vectorize(delta_x)

#define a new method that mixed with numerical and analytical calculation to gain density \xi at z=0
#extend the matter power spectrum
k_pk_arrary = np.logspace(-4,2,1000)
pk_array = vec_power(k_pk_arrary)
split_k = 10 #Mpc^{-1}
for i in range(len(k_pk_arrary)):
    if k_pk_arrary[i] > split_k:
        break
fit_para_rpower, fit_error_rpower = curve_fit(largek_extention, k_pk_arrary[i:-1], pk_array[i:-1])

def xi_dd_z0_mix(r, fit_A, fit_alpha, split_k):
    def integrand(k):
        return power_in_k(k) * k  * np.sin(k*r) / r
    result = integrate.quad(integrand, 0, split_k, epsrel = 1e-5, limit=1000)[0]
    result += (fit_A / (2 * r ** (fit_alpha+3))) * (-1j) ** (-fit_alpha-1)        * (complex(mpmath.gammainc(2 + fit_alpha, -1j * split_k * r)) + (-1) ** (- float(fit_alpha) - 1) * complex(mpmath.gammainc(2 + fit_alpha, 1j * split_k * r)))                 
    return  (1.0 / (2 * np.pi * np.pi) * result).real
vec_xi_dd_z0_mix = np.vectorize(xi_dd_z0_mix)


#points follows gaussian distribution at the s_0 = S(M_max)
def p_delta_s0(delta_s0,s0):
    return 1.0 / ((2 * np.pi * s0) ** 0.5) * np.exp( - delta_s0 ** 2 / (2 * s0))

#calculate the ionization fraction at the decceleration phase
def conditional_BMF(m, s, delta_s0, s0, B0, B1): # BMF = m*dn/dm (from (s0,delta_s0) rather than (0,0))
    deriv = -derivative(sigma, m, dx=1e8) #step of derivative is set to be 1e8 solar mass
    conditional_B0 = B0 + B1 * s0 - delta_s0
    EXP_index = -pow(conditional_B0 + B1 * (s - s0),2.0) / (2.0 * (s - s0))
    BMF = pow(2.0/math.pi,0.5) * rho_bar * deriv * conditional_B0 * (s ** 0.5) / ((s - s0) ** 1.5) * np.exp(EXP_index);
    return BMF;

#calculate the BMF parameters of a certain redshift
def PARA_z(z, M_max, zeta_z_func, T_vir, mu):
    s0 = S(M_max)
    B0 = B_0(z, zeta_z_func, T_vir, mu)
    B1 = B_1(z, zeta_z_func, T_vir, mu)
    return [s0, B0, B1]

#calculate the bubble mass function (dn/dm * m)
#for M_bubble < M_max
def BMF(m, PARA): #the redshift information is in PARA
    [s0, B0, B1] = PARA
    delta_cross = B0 + B1 * s0
    s = S(m)
    def integrand_BMF(delta):
        return p_delta_s0(delta, s0) * conditional_BMF(m, s, delta, s0, B0, B1)
    return integrate.quad(integrand_BMF, -6*(s0**0.5), delta_cross, epsrel = 1e-5)[0]
vec_BMF = np.vectorize(BMF)

# the fraction of the bubbles of M_max in the total mass
def Max_Bubble_fraction(PARA):
    [s0, B0, B1] = PARA
    delta_cross = B0 + B1 * s0
    t_cross = delta_cross / ((2 * s0) ** 0.5)
    return (0.5 - 0.5 * special.erf(t_cross))

#input the redshift array and the neutral fraction history
#output the dxH/dz array
def dxH_dz_cal(z_array, xH_array):
    z_array2 = []
    dxH_dz_array = []
    for i in range(len(z_array) - 1):
        z_array2.append( (z_array[i + 1] + z_array[i]) / 2 )
        dxH_dz_array.append( (xH_array[i + 1] - xH_array[i]) / (z_array[i + 1] - z_array[i]) )
    return[z_array2, dxH_dz_array]

#mass fraction of neutral region
def xH_mass(z, zeta_z_func, T_vir, mu, M_max):
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    ionized_fraction = Max_Bubble_fraction(PARA) + integrate.quad(BMF, (zeta_z_func(z) * TtoM(z, T_vir, mu)), 0.999*M_max, args=(PARA), epsrel = 1e-5)[0] / rho_bar
    return (1 - ionized_fraction)
    
#use BMF_func to gain xH_mass
def xH_mass_interp(z, zeta_z_func, T_vir, mu, M_max, BMF_func):
    ionized_fraction = Max_Bubble_fraction(PARA_z(z, M_max, zeta_z_func, T_vir, mu)) + integrate.quad(BMF_func, (zeta_z_func(z) * TtoM(z, T_vir, mu)), 0.999*M_max, args=(z), epsrel = 1e-5)[0] / rho_bar
    return (1 - ionized_fraction)
vec_xH_mass_interp = np.vectorize(xH_mass_interp)

#volume fraction of neutral region
def xH_vol(z, zeta_z_func, T_vir, mu, M_max):
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    [s0, B0, B1] = PARA 
    def integrand(m): # (m / rho_bar / (1 + delta(z,m)) ) * dn/dm 
        return BMF(m, PARA) / rho_bar / (1 + (B0 + B1 * S(m)) * dicke(z)) 
    delta_cross = B0 + B1 * s0;
    def integrand2(delta): # (m_max / rho_bar / (1 + delta)) * p(delta)
        return p_delta_s0(delta, s0) / (1 + delta * dicke(z))
    xHII_volume = integrate.quad(integrand, (zeta_z_func(z) * TtoM(z,T_vir,mu)), 0.999 * M_max, epsrel = 1e-4)[0] + integrate.quad(integrand2, delta_cross, 6 * (s0 ** 0.5), epsrel = 1e-4)[0]
    xH_volume = 1 - xHII_volume
    return xH_volume
vec_xH_vol = np.vectorize(xH_vol)

#BMF is known, maybe from interpratation
def xH_vol_interp(z, zeta_z_func, T_vir, mu, M_max, BMF_func):
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    [s0, B0, B1] = PARA 
    def integrand(m): # (m / rho_bar / (1 + delta(z,m)) ) * dn/dm 
        return BMF_func(m, z) / rho_bar / (1 + (B0 + B1 * S(m)) * dicke(z)) 
    delta_cross = B0 + B1 * s0;
    def integrand2(delta): # (m_max / rho_bar / (1 + delta)) * p(delta)
        return p_delta_s0(delta, s0) / (1 + delta * dicke(z))
    xHII_volume = integrate.quad(integrand, (zeta_z_func(z) * TtoM(z,T_vir,mu)), 0.999 * M_max, epsrel = 1e-4)[0] + integrate.quad(integrand2, delta_cross, 6 * (s0 ** 0.5), epsrel = 1e-4)[0]
    xH_volume = 1 - xHII_volume
    return xH_volume
vec_xH_vol_interp = np.vectorize(xH_vol_interp)

#define functions to calculate the average temprature of 21cm and CO
def T_21_tilde(z): # in muK
    return 27 * ( ((1+z)/10) * (0.15/OMm/hlittle/hlittle)) ** 0.5 * (OMb * hlittle * hlittle / 0.023) *1e3

#the average CO temperature at redshift z   
def T_CO_bar(z, T_vir, mu): # in muK
    delta_c_z = delta_c / dicke(z)
    M_min = TtoM(z, T_vir, mu)
    f_coll = special.erfc(delta_c_z / (pow(2,0.5) * sigma_z0(M_min)))
    return 59.4 * (1+z) ** 0.5 * f_coll

#the average CO bias at redshift z
def Bias_CO(z, T_vir, mu):
    M_min = TtoM(z, T_vir ,mu)
    growth_factor = dicke(z)
    delta_c_z = delta_c / growth_factor
    sigma_m = sigma_z0(M_min)
    f_coll = special.erfc(delta_c_z / (pow(2,0.5) * sigma_m))
    Gamma_value = special.gamma(1.5) * special.gammaincc(1.5, pow(delta_c_z,2) / (2 * sigma_m * sigma_m))
    b_CO = 1 - growth_factor * growth_factor / delta_c + 1.414 * Gamma_value / (np.pi ** 0.5 * delta_c * f_coll)
    return b_CO

#analytical result for ionization fraction from Mcquinn et al. (2005)
def bar_Q(z, M_max, zeta_z_func, T_vir, mu, PARA):
    [s0, B0, B1] = PARA
    delta_cross = B0 + B1 * s0
    M_min = zeta_z_func(z) * TtoM(z, T_vir, mu)
    s_min = S(M_min)
    s_max = S(M_max)
    s_D = s_min - s_max
    def conditional_Q(delta):
        con_B0 = B0 + B1 * s0 - delta
        return 0.5 * np.exp(- 2 * con_B0 * B1) * special.erfc((con_B0 - B1 * s_D) / (2 * s_D) ** 0.5) \
                + 0.5 * special.erfc((con_B0 + B1 * s_D) / (2 * s_D) ** 0.5)
    def integrand(delta):
        return p_delta_s0(delta, s0) * conditional_Q(delta)
    
    Q = integrate.quad(integrand, - 6 * s0 ** 0.5, delta_cross, epsrel = 1e-4)[0] + Max_Bubble_fraction(PARA)
    return Q

#the anverage (1+delta) for the density of the neutral region
def HIrho_over_rho0(z, zeta_z_func, T_vir, mu, M_max):
    xH_mass_z = 1 - bar_Q(z, M_max, zeta_z_func, T_vir, mu, PARA_z(z, M_max, zeta_z_func, T_vir, mu))
    xH_vol_z = xH_vol(z, zeta_z_func, T_vir, mu, M_max)
    return xH_mass_z / xH_vol_z

#calculate the dxH_dz for bar_Q
def dxH_dz(z, M_max, zeta_z_func, T_vir, mu):
    xH1 = 1 - bar_Q(z - 0.01, M_max, zeta_z_func, T_vir, mu, PARA_z(z - 0.01, M_max, zeta_z_func, T_vir, mu))
    xH2 = 1 - bar_Q(z + 0.01, M_max, zeta_z_func, T_vir, mu, PARA_z(z + 0.01, M_max, zeta_z_func, T_vir, mu))
    return (xH2 - xH1) / 0.02

#the reverse function with known history data(z_dxHdz, dxHdz_ana)
def dxHdz_To_z(dxHdz, M_max, zeta_z_func, T_vir, mu, z_dxHdz, dxHdz_ana):
    if (dxHdz >= max(dxHdz_ana) or dxHdz <= 0):
        print('The value should be an positive value smaller than the maximum dxHI/dz of the history')
        return 0
    z_turn = z_dxHdz[dxHdz_ana.index(max(dxHdz_ana))]
    
    def solution(z_floor, z_top):
        while (abs(z_floor - z_top) > 0.01):
            z_mid = (z_floor + z_top) / 2
            if (dxHdz > dxH_dz(z_mid, M_max, zeta_z_func, T_vir, mu)):
                z_floor = z_mid
            elif (dxHdz == dxH_dz(z_mid, M_max, zeta_z_func, T_vir, mu)):
                z_floor = z_mid
                z_top = z_mid
            else:
                z_top = z_mid
        return (z_floor + z_top) / 2
    z1 = solution(6, z_turn)
    z2 = solution(10, z_turn)
    return [z1, z2] #z1 in the dec phase, z2 in the acc phase

#average bias for ionization
def bar_bias_x(z, M_max, zeta_z_func, T_vir, mu):
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    [s0, B0, B1] = PARA
    M_min = zeta_z_func(z) * TtoM(z, T_vir, mu)
    delta_cross = B0 + B0 * s0
    s_min = S(M_min)
    s_max = S(M_max)
    s_D = s_min - s_max
    def conditional_Q(delta):
        con_B0 = B0 + B1 * s0 - delta
        return 0.5 * np.exp(- 2 * con_B0 * B1) * special.erfc((con_B0 - B1 * s_D) / (2 * s_D) ** 0.5) \
                + 0.5 * special.erfc((con_B0 + B1 * s_D) / (2 * s_D) ** 0.5)
    def integrand(delta):
        dp_ddelta = p_delta_s0(delta, s0) * delta / s0
        return dp_ddelta * conditional_Q(delta)
    
    numerator = integrate.quad(integrand, -6 * s0 ** 0.5, delta_cross, epsrel = 1e-4, limit = 100)[0] \
                + p_delta_s0(delta_cross, s0)
    bar_bias = 1 + numerator / bar_Q(z, M_max, zeta_z_func, T_vir, mu, PARA) / dicke(z)
    return bar_bias

#volume of the bubbles separated d awy with R1 and R2 overlap
def V0_R1R2(R1, R2, r): #R1 >= R2
    if R1 < R2: temp = R1; R1 = R2; R2 = temp
    R_prime_cal = lambda R1,R2,r: (R2**2 - ((r*r + R2*R2 - R1*R1) / (2*r)) ** 2 ) ** 0.5
    V_shadow_cal = lambda R,R_p: 2/3 * np.pi * R**3 - np.pi * R**2 * (R*R - R_p*R_p) ** 0.5 \
                                    + 1/3 * np.pi * (R*R - R_p*R_p) ** 1.5 #bubble radius R
    if R1 <= r - R2:
        #print('R1 = %.3g, R2 = %.3g, d = %.3g. There is no overlapping area.'%(R1, R2, d))
        return 0
    elif R1*R1 <= r*r + R2 * R2:
        R_p = R_prime_cal(R1,R2,r)
        V_shadow_R1 = V_shadow_cal(R1, R_p); V_shadow_R2 = V_shadow_cal(R2, R_p)
        return V_shadow_R1 + V_shadow_R2
    elif R1 < r + R2:
        R_p = R_prime_cal(R1,R2,r)
        V_shadow_R1 = V_shadow_cal(R1, R_p); V_shadow_R2 = V_shadow_cal(R2, R_p)
        return V_shadow_R1 + (4/3 * np.pi * R2**3 -  V_shadow_R2)
    elif R1 > r + R2:
        return 4/3 * np.pi * R2**3
    
#antisymmetric cross correlation at redshift z with x1,x2 sperated r12 away along LoS
#assume that the matter distribution freeze at z, only the mass of the bubble increases
def xi_A_HICO(z, r12, zeta_z_func, HIrho_over_rho0_func, BMF_func_z, M_max, T_vir, mu):
    [z1, z2] = cal_z1_z2(z, r12, 0)
    M_min = zeta_z_func(z) * TtoM(z, T_vir, mu)
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    [s0, B0, B1] = PARA
    
    xi = bar_Q(z, M_max, zeta_z_func, T_vir, mu, PARA)
    xi_z1 = bar_Q(z1, M_max, zeta_z_func, T_vir, mu, PARA_z(z1, M_max, zeta_z_func, T_vir, mu))
    xi_z2 = bar_Q(z2, M_max, zeta_z_func, T_vir, mu, PARA_z(z2, M_max, zeta_z_func, T_vir, mu))
    delta_cross = B0 + B1 * s0
    
    def integrand(m):
        delta_z = (B0 + B1 * S(m)) * dicke(z)
        R0 = MtoR(m / (1 + delta_z)); R1 = MtoR(m / (1 + delta_z) * (xi_z1 / xi))
        V0 = V0_R1R2(R0, R1, r12)
        return (BMF_func_z(m) / m) * V0 * delta_z
    def integrand_max(delta):
        delta_z = delta * dicke(z)
        R0 = MtoR(M_max / (1 + delta_z)); R1 = MtoR(M_max / (1 + delta_z) * (xi_z1 / xi))
        V0 = V0_R1R2(R0, R1, r12)
        return p_delta_s0(delta, s0) * (rho_bar / M_max) * V0 * delta_z
    xi1_delta2 = integrate.quad(integrand, M_min, 0.999 * M_max, epsrel = 1e-4)[0] \
                + integrate.quad(integrand_max, delta_cross, 6 * s0 ** 0.5, epsrel = 1e-4, limit = 200)[0]

    def integrand(m):
        delta_z = (B0 + B1 * S(m)) * dicke(z)
        R0 = MtoR(m / (1 + delta_z)); R2 = MtoR(m / (1 + delta_z) * (xi_z2 / xi))
        V0 = V0_R1R2(R0, R2, r12)
        return (BMF_func_z(m) / m) * V0 * delta_z
    def integrand_max(delta):
        delta_z = delta * dicke(z)
        R0 = MtoR(M_max / (1 + delta_z)); R2 = MtoR(M_max / (1 + delta_z)* (xi_z2 / xi))
        V0 = V0_R1R2(R0, R2, r12)
        return p_delta_s0(delta, s0) * (rho_bar / M_max) * V0 * delta_z
    xi2_delta1 = integrate.quad(integrand, M_min, 0.999 * M_max, epsrel = 1e-4)[0] \
                + integrate.quad(integrand_max, delta_cross, 6 * s0 ** 0.5, epsrel = 1e-4, limit = 200)[0]
    
    xi_HIz1_COz2 = - xi1_delta2 * T_CO_bar(z2, T_vir, mu) * Bias_CO(z2, T_vir, mu) \
                    * T_21_tilde(z1) * HIrho_over_rho0_func(z1)
    xi_HIz2_COz1 = - xi2_delta1 * T_CO_bar(z1, T_vir, mu) * Bias_CO(z1, T_vir, mu) \
                    * T_21_tilde(z2) * HIrho_over_rho0_func(z2)
    return 0.5 * (xi_HIz1_COz2 - xi_HIz2_COz1) # in muK^2
    
#compute the symmetric cross-correlation between 21cm and CO(1-0) line
def xi_S_HICO(z, r12, zeta_z_func, HIrho_over_rho0_func, BMF_func_z, M_max, T_vir, mu):
    M_min = zeta_z_func(z) * TtoM(z, T_vir, mu)
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    [s0, B0, B1] = PARA
    delta_cross = B0 + B1 * s0
    
    def integrand(m):
        delta_z = (B0 + B1 * S(m)) * dicke(z)
        R = MtoR(m / (1 + delta_z)); V0 = V0_R1R2(R, R, r12)
        return (BMF_func_z(m) / m) * V0 * delta_z
    def integrand_max(delta):
        delta_z = delta * dicke(z)
        R = MtoR(M_max / (1 + delta_z)); V0 = V0_R1R2(R, R, r12)
        return p_delta_s0(delta, s0) * (rho_bar / M_max) * V0 * delta_z
    xi_HII_delta = integrate.quad(integrand, M_min, 0.999 * M_max, epsrel = 1e-4)[0] \
                + integrate.quad(integrand_max, delta_cross, 6 * s0 ** 0.5, epsrel = 1e-4, limit = 200)[0]
    
    xi_HI_CO = - xi_HII_delta * T_CO_bar(z, T_vir, mu) * Bias_CO(z, T_vir, mu) * T_21_tilde(z) * HIrho_over_rho0_func(z)
    return xi_HI_CO # in muK^2

#compute the auto cross-correlation of 21cm
def xi_auto_21(z, r12, zeta_z_func, HIrho_over_rho0_func, BMF_func_z, M_max, T_vir, mu):
    M_min = zeta_z_func(z) * TtoM(z, T_vir, mu)
    PARA = PARA_z(z, M_max, zeta_z_func, T_vir, mu)
    [s0, B0, B1] = PARA
    delta_cross = B0 + B1 * s0
    Q = bar_Q(z, M_max, zeta_z_func, T_vir, mu, PARA)
    
    def integrand(m):
        delta_z = (B0 + B1 * S(m)) * dicke(z)
        R = MtoR(m / (1 + delta_z)); V0 = V0_R1R2(R, R, r12)
        return (BMF_func_z(m) / m) * V0
    def integrand_max(delta):
        delta_z = delta * dicke(z)
        R = MtoR(M_max / (1 + delta_z)); V0 = V0_R1R2(R, R, r12)
        return p_delta_s0(delta, s0) * (rho_bar / M_max) * V0
    P1 = integrate.quad(integrand, M_min, 0.999 * M_max, epsrel = 1e-4)[0] \
                + integrate.quad(integrand_max, delta_cross, 6 * s0 ** 0.5, epsrel = 1e-4, limit = 200)[0]
    xi_auto_HII = (1 - Q) * P1
    xi_auto_21 = xi_auto_HII * (T_21_tilde(z) ** 2) * (HIrho_over_rho0_func(z) ** 2)
    return xi_auto_21 # in muK^2
    
#compute the symmetric power spectrum from the cross-correlation function
def Pk_S(k, xi_S_HICO_func, lower_limit, upper_limit):
    integrand = lambda r: xi_S_HICO_func(r) * (r * r) * np.sin(k * r) / (k * r)
    Pk = (hlittle ** 3) * 4 * np.pi * integrate.quad(integrand, lower_limit, upper_limit, epsrel = 1e-3)[0]
    return Pk #in muK^2 h^-3 Mpc^3

#compute the auto power spectrum of CO line
def Pk_auto_CO(z, k, T_vir, mu):
    return (hlittle ** 3) * power_in_k(k) * (dicke(z) ** 2) \
            * (Bias_CO(z, T_vir, mu) ** 2) * (T_CO_bar(z, T_vir, mu) ** 2) #in muK^2 h^-3 Mpc^3
    
#smoothing with points evenly located along LoS in a box of 384 Mpc
def xi_A_HICO_smoothing(z, r12, xi_A_HICO_func, BOX_LEN):
    NUM=20
    r_array = np.linspace(BOX_LEN / NUM, BOX_LEN, NUM)
    xi_A_HICO = 0
    for r in r_array:
        [z1, z2] = cal_z1_z2(z, r, 0)
        xi_A_HICO += xi_A_HICO_func(r12, z1) + xi_A_HICO_func(r12, z2)
    xi_A_HICO += xi_A_HICO_func(r12, z)
    return xi_A_HICO / (2*NUM+1)

#Fourier transformation from antisymmetric cross-correlation to the antisymmetric power spectrum
#under the assumptions that both of the functions above follows the form of dipole
def Pk_A_cal(k, xi_A_HICO_func1d, lower_limit, upper_limit):
    Y10 = (3 / (4 * np.pi)) ** 0.5
    tick1 = time.time()
    def integrand(r,theta):
        integrand = xi_A_HICO_func1d(r) * np.cos(theta) * (r ** 2) * np.sin(theta) * np.sin(k * r * np.cos(theta))
        return integrand;
    result = hlittle**3 * (4 * np.pi) * integrate.dblquad(integrand, 0, np.pi / 2, lower_limit, upper_limit, epsrel = 1e-3)[0]
    result /= Y10
    tick2 = time.time()
    #print(r'k = %3.g (h $Mpc^{-1}$), Pk_A = %6.6g cost time: %3.3g seconds'%(k / antisym_func.hlittle, result, (tick2 - tick1)))
    return result #in muK^2 h^-3 Mpc^3 (the dipole of the antisymmetric power spectrum
    
#the dipole of antisymmetric power spectrum with known fitting parameters
def antisym_fit_curve(kh, A_R, n_R, beta_R, alpha_R):
    return A_R * (kh ** (-n_R)) * np.exp(- beta_R * kh ** alpha_R)
vec_antisym_fit_curve = np.vectorize(antisym_fit_curve)

#fitting the dipole values with power law functions
def fit_nR(kh_array, Pk_A_array, MEASURE):
    if (MEASURE == 1): #power law fitting
        def LS_PowerLaw_factor(k_local, A_R, n_R): #fit the power law between kh (0.14~0.26)
            func_form = A_R * k_local ** (-n_R)
            return func_form
        popt, pocv = curve_fit(LS_PowerLaw_factor, kh_array, Pk_A_array, p0 = [Pk_A_array[-1],2])
    if (MEASURE == 2):
        Pk_A_fit = [] #fitting in the loglog figure, increase the effect of large k part
        for i in range(len(Pk_A_array)):
            Pk_A_fit.append(np.abs(Pk_A_array[i]))
        def LS_PowerLaw_factor(k_local, A_R, n_R): #fit the power law between kh (0.14~0.26)
            func_form = A_R * k_local ** (-n_R)
            return np.log10(func_form)
        popt, pocv = curve_fit(LS_PowerLaw_factor, kh_array, np.log10(Pk_A_fit), p0 = [np.log10(Pk_A_fit[-1]),2])
    return [popt[0], popt[1], pocv[0,0] ** 0.5, pocv[1,1] ** 0.5] #[A_R, n_R, errorbar for A_R, errorbar for n_R]


