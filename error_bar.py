#the method follows Eqn. 8 in Tan's draft
#the variance of dipole term itself of the antisymmetric power spectrum itself

import numpy as np
from scipy import interpolate, integrate  
from scipy.stats import chi2
from scipy.special import erfinv, erfcinv
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
    p = chi2.cdf(signal_noise_ratio_square_sum, free_degree)
    significance = erfinv(p)
    if significance == float('inf'):
        p_comp = integrate.quad(lambda x: chi2.pdf(x, free_degree), signal_noise_ratio_square_sum, np.inf, epsrel = 1e-5)[0]
        significance = erfcinv(p_comp)
    return significance

def reliable_k_range(z, Omega_beam_CO, Omega_patch, delta_nu_CO, delta_D, L_min_21, L_max_21, delta_nu_21):
    """
    Omega_beam_CO : [radian]^2, the instant solid-angle on the sky subtented by one beam of CO detection
    Omega_patch : [radian]^2, the solid-angle on the sky subtented by a box
    delta_nu_CO : MHz, the width of single frequency channel of CO detection
    delta_D : Mpc/h, the length of the observation box along the LoS direction
    L_min_21: meters, the minimum baseline length of 21cm observation
    L_max_21: meters, the maximum baseline length of 21cm observation
    delta_nu_21 : MHz, the width of single frequency channel of 21 detection
    """
    k_per_max_CO = 2.*np.pi / (X(z) * np.sqrt(Omega_beam_CO))
    k_per_min_CO = 2.*np.pi / (X(z) * np.sqrt(Omega_patch))
    k_para_max_CO = 2.*np.pi / (Y(z, nu_CO) * delta_nu_CO)
    k_per_max_21 = 2. * np.pi * L_max_21 / lambda_21 / (1 + z) / X(z);
    k_per_min_21 = 2. * np.pi * L_min_21 / lambda_21 / (1 + z) / X(z);
    k_para_max_21 = 2. * np.pi / (Y(z, nu_21) * delta_nu_21)

    k_per_max = min(k_per_max_CO, k_per_max_21);    k_per_min = max(k_per_min_CO, k_per_min_21)
    k_para_max = min(k_para_max_CO, k_para_max_21); k_para_min = 2.*np.pi / delta_D
    return [k_per_min, k_per_max, k_para_min, k_para_max] # in the unit of h/Mpc

def W(k, mu_0, z, nu_rest, Omega_beam, delta_nu):
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
    sigma_perp = X(z)*np.sqrt(Omega_beam) / (2 * np.pi)
    print('sigma_perp =', sigma_perp)
    sigma_para = Y(z, nu_rest)*delta_nu / (2 * np.pi)
    modification = np.exp(k**2*sigma_perp**2) * integrate.quad(lambda mu: np.exp(k**2*(sigma_para**2-sigma_perp**2)*mu**2), mu_0, 1)[0]
    modification /= (1 - mu_0)
    return modification

def error_CO_auto(k, delta_k, k_limits, P_k, z, Tsys, Omega_patch, delta_D, NUM_PATCH, Omega_beam, delta_nu, Nfeeds, t_int):
    """
    Calculate the error bar on CO auto power spectrum
    -----------
    Parameters:
    -----------
    k : h/Mpc
    delta_k : h/Mpc, the width of the k-shell
    k_limits : h/Mpc, the reliable range [k_per_min, k_per_max, k_para_min, k_para_max]
    P_k : [muK]^2[Mpc/h]^3
    z : redshift
    Tsys : muK
    Omega_patch : [radian]^2, the solid-angle on the sky subtented by a box
    NUM_PATCH : number of patches in a survey
    Omega_beam : [radian]^2, the instant solid-angle on the sky subtented by one beam
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
    # reliable limits in k-space in the unit of h/Mpc
    Omega_survey = NUM_PATCH * Omega_patch
    V_survey = X(z)**2*Omega_survey*delta_D
    t_int *= 1e6
    [k_per_min, k_per_max, k_para_min, k_para_max] = k_limits
    
    if k < (k_para_min**2 + k_per_min**2) ** 0.5:
        print('ERROR: k is below the reliable range')
        return 'ERROR'
    elif k <= k_per_max:
        Nmodes = k**2*delta_k*V_survey/4/np.pi**2
        noise = W(k, 0, z, nu_CO, Omega_beam, delta_nu)*X(z)**2*Y(z, nu_CO)*Tsys**2*Omega_survey/Nfeeds/t_int/np.sqrt(Nmodes)
        cosmic_variance = P_k/np.sqrt(Nmodes)
        total = noise+cosmic_variance
        print('1/W = %4.4g'%W(k, 0, z, nu_CO, Omega_beam, delta_nu))
        print('PN = %4.4g'%(X(z)**2*Y(z, nu_CO)*Tsys**2*Omega_survey/Nfeeds/t_int))
    else:
        mu_0 = (1 - k_per_max**2 / k**2) ** 0.5
        Nmodes = V_survey / (2*np.pi)**3 * 2 * np.pi * k**2 * delta_k * (1 - mu_0)
        noise = W(k, mu_0, z, nu_CO, Omega_beam, delta_nu)*X(z)**2*Y(z, nu_CO)*Tsys**2*Omega_survey/Nfeeds/t_int/np.sqrt(Nmodes)
        cosmic_variance = P_k/np.sqrt(Nmodes)
        total = noise+cosmic_variance
        print(W(k, mu_0, z, nu_CO, Omega_beam, delta_nu))
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

def error_21_auto(k, delta_k, k_limits, P_k, z, Tsys, Omega_patch, delta_D, NUM_PATCH, Omega_beam, Nfeeds, t_int):
    """
    Calculate the error bar on 21cm auto power spectrum
    -----------
    Parameters:
    -----------
    k : h/Mpc
    delta_k : h/Mpc, the width of the k-shell
    k_limits : h/Mpc, the reliable range [k_per_min, k_per_max, k_para_min, k_para_max]
    P_k : [muK]^2[Mpc/h]^3
    z : redshift
    Tsys : muK
    delta_k : h/Mpc, the width of the k-shell 
    Omega_patch : [radian]^2, the solid-angle on the sky subtented by a box
    NUM_PATCH : number of patches in a survey
    Omega_beam : [radian]^2, the instant solid-angle on the sky subtented by one beam
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
    # the k-shell 
    k_inner = k - 0.5 * delta_k
    k_outer = k + 0.5 * delta_k
    # the resolution in k-space
    [k_per_min, k_per_max, k_para_min, k_para_max] = k_limits
    
    delta_k_para = 0.02; delta_k_perp = 0.02 # test
    #delta_k_para = k_para_min
    #delta_k_perp = k_per_min
    
    para_bins_upper = np.int( (k_para_max - delta_k_para) / delta_k_para )
    N_k_para_bin_max = min(np.int(k_outer/delta_k_para - 0.5) + 1, para_bins_upper)
    perp_bins_upper = np.int( (k_per_max - delta_k_perp) / delta_k_perp )
    N_k_perp_bin_max = min(np.int(k_outer/delta_k_perp - 0.5) + 1, perp_bins_upper)
    Omega_survey = NUM_PATCH * Omega_patch
    V_survey = X(z)**2*Omega_survey*delta_D
    t_int *= 1e6 
    total_sum_in_mu, noise_sum_in_mu, cosmic_variance_sum_in_mu, Nmodes = 0, 0, 0, 0
    for i in range(0, N_k_para_bin_max + 1):
            for j in range(0, N_k_perp_bin_max + 1):
                k_para = (i+0.5)*delta_k_para
                k_perp = (j+0.5)*delta_k_perp
                if k_inner <= np.sqrt(k_perp**2+k_para**2) and np.sqrt(k_perp**2+k_para**2) <= k_outer: 
                    u_perp = k_perp*X(z)/2/np.pi
                    Npixel = k_perp * delta_k_perp * delta_k_para * V_survey /4./np.pi**2
                    noise_k_mu = X(z)**2*Y(z, nu_21)*Tsys**2*Omega_patch*NUM_PATCH*Omega_beam/t_int/Nfeeds/n_u(u_perp, lambda_21*(1+z))
                    delta_P_k_mu = P_k + noise_k_mu
                    total_sum_in_mu +=  Npixel/delta_P_k_mu**2
                    cosmic_variance_sum_in_mu +=  Npixel/P_k**2
                    noise_sum_in_mu += Npixel/noise_k_mu**2
                    Nmodes += Npixel

    total = total_sum_in_mu**(-0.5)
    noise = noise_sum_in_mu**(-0.5)
    cosmic_variance = cosmic_variance_sum_in_mu**(-0.5)
    return total, noise, cosmic_variance, Nmodes

def error_21_CO_cross(k, delta_k, k_limits, P_k, z, Omega_patch, NUM_PATCH, delta_D):
    """
    -----------
    Parameters:
    -----------
    k : h/Mpc
    delta_k : h/Mpc, the width of the k-shell
    k_limits : h/Mpc, the reliable range [k_per_min, k_per_max, k_para_min, k_para_max]
    P_k : [muK]^2[Mpc/h]^3
    z : redshift
    Omega_patch : [radian]^2, the solid-angle on the sky subtented by a box
    NUM_PATCH : number of patches in a survey
    delta_D : Mpc/h, the length of the observation box along the LoS direction
    """
    [k_per_min, k_per_max, k_para_min, k_para_max] = k_limits
    V_survey = X(z)**2 * Omega_patch * NUM_PATCH * delta_D
    
    if k < (k_para_min**2 + k_per_min**2) ** 0.5:
        print('ERROR: k is below the reliable range')
        return 'ERROR'
    elif k <= k_per_max:
        Nmodes = k**2*delta_k*V_survey/4/np.pi**2
        error = P_k/np.sqrt(Nmodes)
        #print(k, delta_k, V_survey, Nmodes)
        return error, Nmodes
    else:
        mu_0 = (1 - k_per_max**2 / k**2) ** 0.5
        Nmodes = V_survey / (2*np.pi)**3 * 2 * np.pi * k**2 * delta_k * (1 - mu_0)
        error = P_k/np.sqrt(Nmodes)
        return abs(error), Nmodes

def VAR_PA_k(k, delta_k, z, Omega_patch, NUM_PATCH, Tsys_CO, Tsys_21, Omega_beam_CO, Omega_beam_21, delta_nu_CO, delta_nu_21, Nfeeds_CO, Nfeeds_21, t_int_CO, t_int_21, delta_D, L_min_21, L_max_21, P_k_CO, P_k_21, PS_k, PA_k):
    """
    Calculate the variance on the dipole of the CO auto antisymmetric power spectrum at certain k
    -----------
    Parameters:
    -----------
    k : h/Mpc
    delta_k : h/Mpc, the width of the k-shell
    z : the redshift
    z : redshift
    
    Omega_patch : [radian]^2, the solid-angle on the sky subtented by a box
    NUM_PATCH : number of patches in a survey
    Tsys_CO, Tsys_21 : muK, the system temperature of the two lines detction respectively
    Omega_beam_CO, Omega_beam_21 : radians^2, the solid angle of CO and 21cm detection
    delta_nu_CO, delta_nu_21 : MHz, the width of single frequency channel of CO and 21cm lines
    N_feeds_CO, N_feeds_21 : the number of independent feeds, like different beams, polarizations ...
    t_int_CO, t_int_21 : second, the total integration time
    delta_D : Mpc/h, the length of the observation box along the LoS direction
    L_min_21, L_max_21: meters, the minimum and the maximum baseline length of 21cm observation
    P_k_CO, P_k_21, PS_k, PA_k : [muK]^2[Mpc/h]^3, the auto power spectrum of CO line, the auto power spectrum of 21cm line, the symmetric cross power spectrum and the antisymmetric power spectru,
    --------
    Returns:
    --------
    sigma_error: [muK]^2[Mpc/h]^3, the total error the dipole detection
    sigma_cosmic_variance : [muK]^2[Mpc/h]^3, the cosmic variance of the dipole detection
    """

    k_limits = reliable_k_range(z, Omega_beam_CO, Omega_patch, delta_nu_CO, delta_D, L_min_21, L_max_21, delta_nu_21)
    k_perp_min, k_perp_max, k_para_min, k_para_max = k_limits
    if k**2 < k_perp_min**2 + k_para_min**2 or k**2 > k_perp_max**2 + k_para_max**2:
        return 'ERROR: k is out of the reliable range'
    
    sigma_P_k_CO, _, sigma_P_k_CO_cv, _ = error_CO_auto(k, delta_k, k_limits, P_k_CO, z, Tsys_CO, Omega_patch, delta_D, NUM_PATCH, Omega_beam_CO, delta_nu_CO, Nfeeds_CO, t_int_CO)
    sigma_P_k_21, _, sigma_P_k_21_cv, _ = error_21_auto(k, delta_k, k_limits, P_k_21, z, Tsys_21, Omega_patch, delta_D, NUM_PATCH, Omega_beam_21, Nfeeds_21, t_int_21)
    sigma_PS_k, _ = error_21_CO_cross(k, delta_k, k_limits, PS_k, z, Omega_patch, NUM_PATCH, delta_D)
    sigma_PA_k, _ = error_21_CO_cross(k, delta_k, k_limits, PA_k, z, Omega_patch, NUM_PATCH, delta_D) #cosmic variance only
    VAR_PA_k = (9 / 20 / np.pi) * sigma_PA_k**2 + sigma_P_k_21 * sigma_P_k_CO - sigma_PS_k**2
    VAR_PA_k_cv = (9 / 20 / np.pi) * sigma_PA_k**2 + sigma_P_k_21_cv * sigma_P_k_CO_cv - sigma_PS_k**2 #eqn 9 in Tan's draft
    print('sigma_P_k_CO = %5.5g'%sigma_P_k_CO,'sigma_P_k_CO_cv = %5.5g'%sigma_P_k_CO_cv)
    print('sigma_P_k_21 = %5.5g'%sigma_P_k_21,'sigma_P_k_21_cv = %5.5g'%sigma_P_k_21_cv)
    print('sigma_PS_k = %5.5g'%sigma_PS_k,'sigma_PA_k = %5.5g'%sigma_PS_k)
    print('sigma_error = %5.5g'%VAR_PA_k**0.5,'sigma_cv = %5.5g'%VAR_PA_k_cv**0.5)
    return VAR_PA_k**0.5, VAR_PA_k_cv**0.5 #return the sigma of error
