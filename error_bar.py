import numpy as np
from scipy import interpolate, integrate  
from scipy.stats import chi2
from scipy.special import erfinv
#*****************************#
## Cosmological Parameters ##
#*****************************#
HubblePara = 0.678
## no unit
HubbleLength = 3000
## in unit Mpc/h

OMega_m = 0.306
OMega_r = 8.6e-5
OMega_lambda = 1.0-OMega_m-OMega_r
## Constants of LambdaCDM, no units

nu_21 = 1420.405752
nu_CO = 115270
## Frequency (MHz) of emission lines at rest

lambda_21 = 0.2110611405
lambda_CO = 0.00261
## Wavelength (meters) of emission lines at rest

def E(z):
    """
    The scale factor which relates the current Hubble parameter H0 with the Hubble parameter at redshift z H(z) by H(z) = E(z) H0. 
    """
    return np.sqrt( OMega_m*(1+z)**3 + OMega_r*(1+z)**4 + OMega_lambda )

def X(z):
    """
    The radial co-moving distance (Mpc/h) at redshift z.
    """
    integral = integrate.quad(lambda x: 1./E(x), 0, z)[0]
    integral *= HubbleLength
    return integral

def Y(z, nu_rest):
    """
    Y = d r_parallel (Mpc/h) / d v (MHz) = Hubblelength * (1+z)^2/E(z)/nu_rest.
    
    """
    return HubbleLength * (1.+z) **2. / E(z) / nu_rest 

#*****************************#
## Main Functions ##
#*****************************#

def significance_level(signal_noise_ratio_square_sum, free_degree):
    """
    Calculate the significance level for claiming a detection.
    Assume the square sum of signal-noise-ratio from multiple measurements follows a chi-square distribution,
    then get the p-value and convert it to the equivalent significance level of a normal distribution. 
    """
    p = chi2.cdf(square_sum, free_degree)
    return erfinv(p)

def W(k, z, nu_rest, Omega_beam, delta_nu):
    """
    The modified factor due to limited resolution
    -----------
    Parameters:
    -----------
    k : h/Mpc
    z : redshift
    nu_rest : MHz, the frequency of the emission line at its rest frame
    Omega_beam : [radian]^2, the instant solid-angle on the sky subtented by one beam 
    delta_nu : MHz, the width of single frequency channel
    """
    sigma_perp = X(z)*np.sqrt(Omega_beam)
    sigma_para = Y(z, nu_rest)*delta_nu
    modification = np.exp(k**2*sigma_perp**2)/integrate.quad(lambda mu: np.exp(-k**2*(sigma_para**2-sigma_perp**2)*mu**2), 0, 1)[0]
    return modification

def error_CO_auto(k, P_k, z, Tsys, delta_k, Omega_survey, Omega_beam, B, delta_nu, Nfeeds, t_int):
    """
    Calculate the error bar on CO auto power spectrum
    -----------
    Parameters:
    -----------
    k : h/Mpc
    P_k : [muK]^2[Mpc/h]^3
    z : redshift
    Tsys : muK
    delta_k : h/Mpc, the width of the k-shell 
    Omega_survey : [radian]^2, the solid-angle on the sky subtented by the whole survey
    Omega_beam : [radian]^2, the instant solid-angle on the sky subtented by one beam 
    B : MHz, the full bandwidth of the receiver
    delta_nu : MHz, the width of single frequency channel
    Nfeeds : the number of independent feeds, like different beams, polarizations ...
    t_int : second, the total integration time
    --------
    Returns:
    --------
    total : [muK]^2[Mpc/h]^3
    noise : [muK]^2[Mpc/h]^3
    cosmic_variance : [muK]^2[Mpc/h]^3
    Nmodes : 
    """
    V_survey = X(z)**2*Y(z, nu_21)*Omega_survey*B
    t_int *= 1e6
    Nmodes = k**2*delta_k*V_survey/4/np.pi**2
    noise = W(k, z, nu_CO, Omega_beam, delta_nu)*X(z)**2*Y(z, nu_CO)*Tsys**2*Omega_survey/Nfeeds/t_int/np.sqrt(Nmodes)
    cosmic_variance = W(k, z, nu_CO, Omega_beam, delta_nu)*P_k/np.sqrt(Nmodes)
    total = noise+cosmic_variance

    return total, noise, cosmic_variance, Nmodes

def n_u(u, wavelength):
    """
    The number density of redundant baselines in uv space, assuming a SKA-low core configuration.
    u : no unit, L/wavelength
    wavelength : m
    """
    baseline_number_density = [0.0609823,0.05809261,0.05520293,0.05231325,0.04942356,\
    0.04673317,0.04404277,0.04135238,0.03866198,0.03597158,0.03328119,0.03079008,\
    0.02829897,0.02580787,0.02331676,0.02102494,0.01873312,0.01644131,0.01414949,\
    0.01205696,0.00996443,0.00807119,0.00617795,0.00468328,0.00318862,0.00209253,0.00099644,0.00049822,0.]
    baseline_length = [ 35, 70, 105,  140,  175,  210,  245,  280,  315,  350,  385, 420,  455,  \
    490,  525,  560,  595,  630,  665,  700,  735,  770, 805,  840,  875,  910,  945,  980, 1015 ]
    f = interpolate.interp1d(baseline_length, baseline_number_density, kind = "cubic", fill_value="extrapolate")
    return f(u*wavelength)*wavelength**2

def error_21_auto(k, P_k, z, Tsys, delta_k, Omega_survey, Omega_beam, B, delta_nu, Nfeeds, t_int, Lmin):
    """
    Calculate the error bar on 21cm auto power spectrum
    -----------
    Parameters:
    -----------
    k : h/Mpc
    P_k : [muK]^2[Mpc/h]^3
    z : redshift
    Tsys : muK
    delta_k : h/Mpc, the width of the k-shell 
    Omega_survey : [radian]^2, the solid-angle on the sky subtented by the whole survey
    Omega_beam : [radian]^2, the instant solid-angle on the sky subtented by one beam 
    B : MHz, the full bandwidth of the receiver
    delta_nu : MHz, the width of single frequency channel
    Nfeeds : the number of independent feeds, like different beams, polarizations ...
    t_int : second, the total integration time
    Lmin : m, minimum baseline length
    --------
    Returns:
    --------
    total : [muK]^2[Mpc/h]^3
    noise : [muK]^2[Mpc/h]^3
    cosmic_variance : [muK]^2[Mpc/h]^3
    Nmodes : 
    """
    # the k-shell 
    k_inner = k
    k_outer = k+delta_k
    # the resolution in k-space
    delta_k_para = 2*np.pi/Y(z, nu_21)/B
    delta_k_perp = 2*np.pi*Lmin/lambda_21/(1.+z)/X(z)

    N_k_para_bin_max = np.int(k_outer/delta_k_para - 0.5) + 1
    N_k_perp_bin_max = np.int(k_outer/delta_k_perp - 0.5) + 1
    V_survey = X(z)**2*Y(z, nu_21)*Omega_survey*B
    t_int *= 1e6 
    total_sum_in_mu, noise_sum_in_mu, cosmic_variance_sum_in_mu, Nmodes = 0, 0, 0, 0
    for i in range(0, N_k_para_bin_max + 1):
            for j in range(0, N_k_perp_bin_max + 1):
                k_para = (i+0.5)*delta_k_para
                k_perp = (j+0.5)*delta_k_perp
                if k_inner <= np.sqrt(k_perp**2+k_para**2) and np.sqrt(k_perp**2+k_para**2) <= k_outer: 
                    u_perp = k_perp*X(z)/2/np.pi
                    Npixel = k_perp * delta_k_perp * delta_k_para * V_survey /4./np.pi**2
                    noise_k_mu = X(z)**2*Y(z, nu_21)*Tsys**2*Omega_survey*Omega_beam/t_int/Nfeeds/n_u(u_perp, lambda_21*(1+z))
                    delta_P_k_mu = P_k + noise_k_mu
                    total_sum_in_mu +=  Npixel/delta_P_k_mu**2
                    cosmic_variance_sum_in_mu +=  Npixel/P_k**2
                    noise_sum_in_mu += Npixel/noise_k_mu**2
                    Nmodes += Npixel

    total = total_sum_in_mu**(-0.5)
    noise = noise_sum_in_mu**(-0.5)
    cosmic_variance = cosmic_variance_sum_in_mu**(-0.5)
    return total, noise, cosmic_variance, Nmodes

def error_21_CO_cross(k, P_k, z, delta_k, Omega_survey, B):
    """
    -----------
    Parameters:
    -----------
    k : h/Mpc
    P_k : [muK]^2[Mpc/h]^3
    z : redshift
    delta_k : h/Mpc, the width of the k-shell 
    Omega_survey : [radian]^2, the solid-angle on the sky subtented by the whole survey
    B : MHz, the full bandwidth of the receiver
    """
    V_survey = X(z)**2*Y(z, nu_21)*Omega_survey*B
    N_modes = V_survey*k**2*delta_k/4/np.pi**2
    error = P_k/np.sqrt(N_modes) 

    return error, N_modes          







