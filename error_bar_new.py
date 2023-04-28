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

T_sys_CO = 25e6
## system tempertature (muK) of CO observations

def E(z):
    """
    The scale factor which relates the current Hubble parameter H0 with the Hubble parameter at redshift z H(z) by H(z) = E(z) H0. 
    """
    return np.sqrt( OMega_m*(1+z)**3 + OMega_r*(1+z)**4 + OMega_lambda )

def X(z):
    """
    The radial co-moving distance (Mpc/h) at redshift z.
    ---------------------
    INPUT:
    z : the redshift we are interested
    ---------------------
    OUTPUT:
    X : Mpc/h, X * theta ~ r_perp (Mpc/h)
    """
    integral = integrate.quad(lambda x: 1./E(x), 0, z)[0]
    integral *= HubbleLength
    return integral

def Y(z, nu_rest):
    """
    Y = d r_parallel (Mpc/h) / d v (MHz) = Hubblelength * (1+z)^2/E(z)/nu_rest.
    ---------------------
    INPUT:
    z : the redshift we are interested
    nu_rest : Mhz, the rest-frame frequency of the line we are interested in
    ---------------------
    OUTPUT:
    Y : Mpc/h/Mhz, Y * dv(MHz) ~ r_parallel(Mpz/h)
    """
    return HubbleLength * (1.+z) **2. / E(z) / nu_rest
    
#***************************************#
## Telescope Parameter Computation ##
#***************************************#

def core_diameter(num_station, station_diameter):
    '''
    compute the core size of the SKA2-low, assuming stations following the densest packing--hexagonal packing
    ---------------------
    INPUT:
    num_circle : the number of circles (the number of stations in SKA2-low)
    ---------------------
    OUTPUT:
    core_diameter : the diameter of core area, in which all of the stations can be put in
    '''
    eta = 0.9069 # the hexagonal packing density
    station_area = num_station * (np.pi * station_diameter * station_diameter / 4)
    core_area = station_area / eta
    core_diameter = 2 * (core_area / np.pi)**0.5
    return core_diameter

def baseline_computation(d_station, N_stations, r_core):
    '''
    compute the number density of SKA-low stations along one certain direction
    ---------------------
    INPUT:
    d_station : the diamter of the stations
    N_stations : the number of stations
    r_core : m,  the radius of the SKA-low core
    ---------------------
    OUTPUT:
    baseline_length : m, the baseline length
    baseline_number_density : the number density corresponding to the baseline along one direcntion
    '''
    station_density = N_stations/np.pi/r_core**2
    Nsides = np.int(2*r_core/d_station)+1
    weights_pixel = np.zeros((Nsides+1, Nsides+1))
    for i in range(0, Nsides+1):
        for j in range(0, Nsides+1):
            x_cor, y_cor = (i+0.5), (j+0.5)
            if np.sqrt((x_cor-Nsides/2)**2+(y_cor-Nsides/2)**2)*d_station <= r_core:
                weights_pixel[i,j] = 1
    baseline_length, baseline_number_density = np.arange(1, Nsides+1), np.zeros(Nsides)
    for l in baseline_length:
        for i in range(0, Nsides+1):
            for j in range(0, Nsides+1):
                if( j+l <= Nsides ):
                    baseline_number_density[l-1] += weights_pixel[i,j]*weights_pixel[i,j+l]*station_density**2*d_station**2
    baseline_length *= d_station
    #baseline_number_density*np.pi*r_core**2
    return baseline_length, baseline_number_density

def FoV_computation(z, lambda_rest, D):
    '''
    compute the antenna angular resolution
    ---------------------
    INPUT:
    z : the redshift we are interested
    lambda_rest : m, the restframe wavelength
    D : n, the diamter of the antenna
    ---------------------
    OUTPUT:
    Omega_FoV : radians^2, antenna field of view area (single beam FoV)
    '''
    return np.pi / 4 * (1.3 * lambda_rest * (1+z) / D)**2

def angular_resolution(z, lambda_rest, D):
    '''
    compute the antenna angular resolution
    ---------------------
    INPUT:
    z : the redshift we are interested
    lambda_rest : m, the restframe wavelength
    D : m, the diamter of the antenna
    ---------------------
    OUTPUT:
    theta_min : radians, the angular resolution
    '''
    return 1.22 * lambda_rest * (1 + z) / D

def effective_area_21(z, D):
    '''
    compute the effective area of a SKA-low station
    ---------------------
    INPUT:
    z : the redshift we are interested
    D : m, the diamter of the antenna
    ---------------------
    OUTPUT:
    A_e : m^2, the effective area of a SKA-low station
    '''
    nu = 1420.405752 / (1 + z) #Mhz
    A_e = np.pi * (D / 2)**2 * (94 / nu)**2
    return A_e

def beam_21_computation(z, D):
    '''
    compute the beam size of a single SKA-low station
    --------------------
    INPUT:
    z : the redshift we are interested
    D : m, the diamter of the antenna
    ---------------------
    OUTPUT:
    Omega_beam : radians^2, antenna beam area (single station)
    '''
    return (lambda_21 * (1 + z))**2 / effective_area_21(z, D)

def beam_CO_computation(z, A_eff):
    '''
    compute the beam size of the SKA-mid dish
    --------------------
    INPUT:
    z : the redshift we are interested
    A_eff : m^2, the effective area of a SKA-mid dish
    ---------------------
    OUTPUT:
    Omega_beam : radians^2, antenna beam area (single dish)
    '''
    return (lambda_CO * (1 + z))**2 / A_eff

def bandwidth_computation(z, nu_rest, BOX_LEN):
    '''
    compute the bandwidth of the correspondent survey LoS length at redshift z
    ---------------------
    INPUT:
    z : the redshift we are interested
    BOX_LEN : Mpc/h, the survey depth along the LoS direction (the length of the box side)
    ---------------------
    OUTPUT:
    B : Mhz, the bandwidth of the observation
    '''
    return BOX_LEN / Y(z, nu_rest)
    
def sky_area_patch(z, BOX_LEN):
    '''
    compute the sky area covered by a single observation box of side length BOX_LEN
    ---------------------
    INPUT:
    z : the redshift we are interested
    BOX_LEN : Mpc/h, the survey depth along the LoS direction (the length of the box side)
    ---------------------
    OUTPUT:
    Omega_patch : [radians]^2, the sky area covered by a single observation box
    '''
    return BOX_LEN**2 / X(z)**2

def T_sys_21(z):
    '''
    compute the system temperature of EoR 21cm observations
    ---------------------
    INPUT:
    z : the redshift we are interested
    ---------------------
    OUTPUT:
    T_sys_21 : muK, the system temperature of EoR 21cm observations
    '''
    return 280 * ((1 + z) / 7.4)**2.3 * 1e6
#***************************************#
## Main Functions ##
#***************************************#

def reliable_k_range(z, D_dish, BOX_LEN, NUM_PATCH, delta_nu_CO, L_min_21, L_max_21, delta_nu_21):
    """
    The modified factor due to limited resolution
    -----------
    Parameters:
    -----------
    z : redshift we are interested in
    D_dish : m, the diameter of the SKA-mid dish
    BOX_LEN : Mpc/h, the box length of a single observation box
    NUM_PATCH : the number of boxes
    delta_nu_CO : MHz, the width of single frequency channel of CO detection
    L_min_21: meters, the minimum baseline length of 21cm observation
    L_max_21: meters, the maximum baseline length of 21cm observation
    delta_nu_21 : MHz, the width of single frequency channel of 21 detection
    --------
    Returns:
    --------
    [k_per_min, k_per_max, k_para_min, k_para_max] : h/Mpc, the limit of the reliable
        k modes
    """
    theta_min_CO = angular_resolution(z, lambda_CO, D_dish)
    Omega_survey = sky_area_patch(z, BOX_LEN) * NUM_PATCH
    
    k_per_max_CO = 2.*np.pi / (X(z) * theta_min_CO)
    k_per_min_CO = 2.*np.pi / (X(z) * np.sqrt(Omega_survey))
    k_para_max_CO = 2.*np.pi / (Y(z, nu_CO) * delta_nu_CO)
    k_per_max_21 = 2. * np.pi * L_max_21 / lambda_21 / (1 + z) / X(z);
    k_per_min_21 = 2. * np.pi * L_min_21 / lambda_21 / (1 + z) / X(z);
    k_para_max_21 = 2. * np.pi / (Y(z, nu_21) * delta_nu_21)

    k_per_max = min(k_per_max_CO, k_per_max_21);    k_per_min = max(k_per_min_CO, k_per_min_21)
    k_para_max = min(k_para_max_CO, k_para_max_21); k_para_min = 2.*np.pi / BOX_LEN
    return [k_per_min, k_per_max, k_para_min, k_para_max] # in the unit of h/Mpc

def sphere_k_modes(k, delta_k, k_limits, BOX_LEN, NUM_PATCH):
    """
    Compute the number of k modes in the sphere within the limits of reliable k mode
    -----------
    Parameters:
    -----------
    k : h/Mpc
    Delta_k : h/Mpc, the width of the k-shell
    k_limits : [k_per_min, k_per_max, k_para_min, k_para_max], h/Mpc, the limit of
        the reliable k modes
    BOX_LEN : Mpc/h, the box length of a single observation box
    NUM_PATCH : the number of boxes
    -----------
    Return:
    -----------
    N_modes : the number of reliable k modes in the limits
    """
    k_per_min, k_per_max, k_para_min, k_para_max = k_limits
    mu1 = max(k**2 - k_per_max**2, 0)**0.5 / k;  mu2 = max(k**2 - k_per_min**2, 0)**0.5 / k
    mu3 = min(k_para_min / k, 1);                mu4 = min(k_para_max / k, 1)
    
    V_survey = BOX_LEN**3 * NUM_PATCH
    n_modes = V_survey / (2 * np.pi)**3
    if mu2 <= mu3:
        print('ALLERT : usually this case not happening')
        N_modes = n_modes * 2 * np.pi * k**2 * delta_k * (mu4 - mu3 + mu2 - mu1)
    else:
        N_modes = n_modes * 2 * np.pi * k**2 * delta_k * (mu4 - mu1)
    return N_modes

def W_CO(z, k, k_limits, D_dish, delta_nu_CO):
    """
    The modified factor due to limited resolution of SKA-mid dishes
    -----------
    Parameters:
    -----------
    z : redshift
    k : h/Mpc
    k_limits : [k_per_min, k_per_max, k_para_min, k_para_max], h/Mpc, the limit of
        the reliable k modes
    D_dish : m, the diameter of the SKA-mid dish
    delta_nu_CO : MHz, the width of single frequency channel of SKA-mid antenna
    -----------
    Return:
    -----------
    W(k) : the modification factor
    """
    k_per_min, k_per_max, k_para_min, k_para_max = k_limits
    mu1 = max(k**2 - k_per_max**2, 0)**0.5 / k;  mu2 = max(k**2 - k_per_min**2, 0)**0.5 / k
    mu3 = min(k_para_min / k, 1);                mu4 = min(k_para_max / k, 1)
    
    theta_FWHM = 1.22 * (lambda_CO * (1+z) / D_dish)
    sigma_perp = X(z) * theta_FWHM / (8 * np.log(2))**0.5
    print('sigma_perp =', sigma_perp)
    sigma_para = Y(z, nu_CO) * delta_nu_CO #/ 2 / np.pi
    if mu2 <= mu3:
        print('ALLERT : usually this case not happening')
        modification = integrate.quad(lambda mu: np.exp(-k**2 * (sigma_para**2-sigma_perp**2) * mu**2), mu1, mu2, epsrel = 1e-4)[0] + integrate.quad(lambda mu: np.exp(-k**2 * (sigma_para**2-sigma_perp**2) * mu**2), mu3, mu4)[0]
        modification *= np.exp(-k**2 * sigma_perp**2) / (mu4 - mu3 + mu2 - mu1)
    else:
        modification = integrate.quad(lambda mu: np.exp(-k**2 * (sigma_para**2-sigma_perp**2) * mu**2), mu1, mu4, epsrel = 1e-4)[0]
        modification *= np.exp(-k**2 * sigma_perp**2) / (mu4 - mu1)
    return modification

def error_CO_auto(z, k, delta_k, k_limits, P_k, BOX_LEN, NUM_PATCH, D_dish, delta_nu_CO, N_antenna_mid, t_total):
    """
    Calculate the error bar on CO auto power spectrum
    -----------
    Parameters:
    -----------
    z : redshift
    k : h/Mpc
    delta_k : h/Mpc, the width of the k-shell
    k_limits : h/Mpc, [k_per_min, k_per_max, k_para_min, k_para_max], the limits of reliable k range
    P_k : [muK]^2[Mpc/h]^3, the auto CO power spectrum
    BOX_LEN : Mpc/h, the box length of a single observation box
    NUM_PATCH : number of patches in a survey
    D_dish : m, the diameter of the SKA-mid dish
    delta_nu_CO : MHz, the width of single frequency channel of SKA-mid
    N_antenna_mid : the number of SKA-mid antennas
    t_total : hours, the total integration time
    --------
    Returns:
    --------
    total : [muK]^2[Mpc/h]^3
    noise : [muK]^2[Mpc/h]^3
    cosmic_variance : [muK]^2[Mpc/h]^3
    Nmodes : the number of k modes within the shell (k +- delta_k)
    """
    # reliable limits in k-space in the unit of h/Mpc
    
    Omega_survey = NUM_PATCH * sky_area_patch(z, BOX_LEN)
    V_survey = BOX_LEN**3 * NUM_PATCH
    t_total *= 3600 * 1e6 #transform the unit of hour to '1e-6 second'
    N_feeds = N_antenna_mid * 2 # 2 polarizations
    
    P_N = X(z)**2 * Y(z, nu_CO) * T_sys_CO**2 * Omega_survey / (N_feeds * t_total)
    W_k = W_CO(z, k, k_limits, D_dish, delta_nu_CO)
    print('1/W = %4.4g'%(1/W_k))
    print('PN = %4.4g'%P_N)
    
    N_modes = sphere_k_modes(k, delta_k, k_limits, BOX_LEN, NUM_PATCH)
    sigma_cosmic_variance = P_k / (N_modes)**0.5
    sigma_noise = P_N / N_modes**0.5 / W_k
    sigma_error = sigma_cosmic_variance / W_k + sigma_noise
    return sigma_error, sigma_noise, sigma_cosmic_variance, N_modes

def n_u(z, d_station, N_stations, r_core):
    """
    The number density of redundant baselines in uv space, assuming a SKA-low core configuration.
    -----------
    Parameters:
    -----------
    u : no unit, L/wavelength
    z : the redshift we are interested in
    d_station : m, the diamter of the stations
    N_stations : the number of stations
    r_core : m,  the radius of the SKA-low core
    --------
    Returns:
    --------
    n_u_function : n(u), the function of baseline number density in uv plane
    """
    wavelength = lambda_21 * (1 + z)
    baseline_length, baseline_number_density = baseline_computation(d_station, N_stations, r_core)
    f = interpolate.interp1d(baseline_length, baseline_number_density, kind = "cubic", fill_value="extrapolate")
    n_u_function = lambda u: f(u*wavelength)*wavelength**2
    return n_u_function
    
def error_21_auto(z, k, delta_k, k_limits, P_k, BOX_LEN, NUM_PATCH, d_station, N_stations, L_max_21, t_total):
    '''
    calculate the error bar on 21cm auto power spectrum
    -----------
    Parameters:
    -----------
    z : redshift
    k : h/Mpc
    delta_k : h/Mpc, the width of the k-shell
    k_limits : h/Mpc, [k_per_min, k_per_max, k_para_min, k_para_max], the limits of reliable k range
    P_k : [muK]^2[Mpc/h]^3, the auto 21 power spectrum
    BOX_LEN : Mpc/h, the box length of a single observation box
    NUM_PATCH : number of patches in a survey
    d_station : m, the diameter of the SKA-low station
    N_stations : the number of SKA-low stations in the core
    L_max_21: meters, the maximum baseline length of 21cm observation, the diameter of the core
    t_total : hours, the total integration time
    '''
    k_per_min, k_per_max, k_para_min, k_para_max = k_limits
    mu1 = max(k**2 - k_per_max**2, 0)**0.5 / k;  mu2 = max(k**2 - k_per_min**2, 0)**0.5 / k
    mu3 = min(k_para_min / k, 1);                mu4 = min(k_para_max / k, 1)
    
    t_total *= 3600 * 1e6 #transform the unit of hour to '1e-6 second'
    V_survey = BOX_LEN**3 * NUM_PATCH
    Omega_survey = NUM_PATCH * sky_area_patch(z, BOX_LEN)
    Omega_beam = beam_21_computation(z, d_station)
    u_perp_func = lambda mu: X(z) * k / (2 * np.pi) * (1 - mu**2)**0.5
    n_u_func = n_u(z, d_station, N_stations, L_max_21 / 2)
    PN_mu = lambda mu: X(z)**2 * Y(z, nu_21) * (Omega_survey * Omega_beam) \
            / (2 * n_u_func(u_perp_func(mu)) * t_total) * T_sys_21(z)**2  #N_polarizatin = 2
            
    integrand_noise = lambda mu: k**2 * delta_k * V_survey / (4 * np.pi**2 * PN_mu(mu)**2 )
    integrand_total = lambda mu: k**2 * delta_k * V_survey / (4 * np.pi**2 * (PN_mu(mu) + P_k)**2 )
    if mu2 <= mu3:
        print('ALLERT : usually this case not happening')
        N_modes = k**2 * delta_k * V_survey / (4 * np.pi**2) * (mu4 - mu3 + mu2 - mu1)
        sigma_cosmic_variance = P_k / N_modes**0.5
        sigma_noise = (integrate.quad(integrand_noise, mu1, mu2, epsrel = 1e-4)[0] \
                        + integrate.quad(integrand_noise, mu3, mu4, epsrel = 1e-4)[0]) ** (-0.5)
        sigma_total = (integrate.quad(integrand_total, mu1, mu2, epsrel = 1e-4)[0] \
                        + integrate.quad(integrand_total, mu3, mu4, epsrel = 1e-4)[0]) ** (-0.5)
    else:
        N_modes = k**2 * delta_k * V_survey / (4 * np.pi**2) * (mu4 - mu1)
        sigma_cosmic_variance = P_k / N_modes**0.5
        sigma_noise = integrate.quad(integrand_noise, mu1, mu4, epsrel = 1e-4)[0] ** (-0.5)
        sigma_total = integrate.quad(integrand_total, mu1, mu4, epsrel = 1e-4)[0] ** (-0.5)
    #print('21cm auto:',sigma_total, sigma_noise, sigma_cosmic_variance, N_modes)
    return sigma_total, sigma_noise, sigma_cosmic_variance, N_modes
    
def error_21_CO_cross(k, delta_k, k_limits, P_k, BOX_LEN, NUM_PATCH):
    """
    compute the error bar of the symmetric and anti-symmetric cross power spectrum
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
    --------
    Returns:
    --------
    total : [muK]^2[Mpc/h]^3
    noise : [muK]^2[Mpc/h]^3
    cosmic_variance : [muK]^2[Mpc/h]^3
    Nmodes : the number of k modes within the shell (k +- delta_k)
    """
    N_modes = sphere_k_modes(k, delta_k, k_limits, BOX_LEN, NUM_PATCH)
    sigma_error = abs(P_k) / (N_modes)**0.5
    return sigma_error, N_modes

def error_PA_k(z, k, delta_k, t_total_CO, t_total_21, BOX_LEN, NUM_PATCH, delta_nu_CO, delta_nu_21, D_dish, N_antenna_mid, d_station, N_stations, L_min_21, L_max_21, P_k_CO, P_k_21, PS_k, PA_k):
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
    t_int_CO, t_int_21 : hour, the total integration time of CO and 21cm observation respectively
    BOX_LEN : Mpc/h, the side length of a observation box
    NUM_PATCH : the number of observation boxes
    delta_nu_CO, delta_nu_21 : MHz, the width of single frequency channel of CO and 21cm lines
    D_dish : m, the diameter of the SKA-mid dish
    N_antenna_mid : the number of SKA-mid antennas
    d_station : m, the diameter of the SKA-low station
    N_stations : the number of SKA-low stations in the core
    L_min_21, L_max_21: meters, the minimum and the maximum baseline length of 21cm observation
    P_k_CO, P_k_21, PS_k, PA_k : [muK]^2[Mpc/h]^3, the auto power spectrum of CO line, the auto power spectrum of 21cm line, the symmetric cross power spectrum and the antisymmetric power spectru,
    --------
    Returns:
    --------
    sigma_error: [muK]^2[Mpc/h]^3, the total error the dipole detection
    sigma_cosmic_variance : [muK]^2[Mpc/h]^3, sigma of the cosmic variance of the dipole detection
    """

    k_limits = reliable_k_range(z, D_dish, BOX_LEN, NUM_PATCH, delta_nu_CO, L_min_21, L_max_21, delta_nu_21)
    k_per_min, k_per_max, k_para_min, k_para_max = k_limits
    if k**2 < k_per_min**2 + k_para_min**2 or k**2 > k_per_max**2 + k_para_max**2:
        return 'ERROR: k is out of the reliable range'
    
    sigma_P_k_CO, _, sigma_P_k_CO_cv, N_modes = error_CO_auto(z, k, delta_k, k_limits, P_k_CO, BOX_LEN, NUM_PATCH, D_dish, delta_nu_CO, N_antenna_mid, t_total_CO)
    sigma_P_k_21, _, sigma_P_k_21_cv, N_modes = error_21_auto(z, k, delta_k, k_limits, P_k_21, BOX_LEN, NUM_PATCH, d_station, N_stations, L_max_21, t_total_21)
    sigma_PS_k, _ = error_21_CO_cross(k, delta_k, k_limits, PS_k, BOX_LEN, NUM_PATCH)
    sigma_PA_k, _ = error_21_CO_cross(k, delta_k, k_limits, PA_k, BOX_LEN, NUM_PATCH)
    
    print('sigma_P_k_CO = %5.5g'%sigma_P_k_CO, 'sigma_P_k_CO_cv = %5.5g'%sigma_P_k_CO_cv)
    print('sigma_P_k_21 = %5.5g'%sigma_P_k_21, 'sigma_P_k_21_cv = %5.5g'%sigma_P_k_21_cv)
    print('sigma_PS_k = %5.5g'%sigma_PS_k, 'sigma_PA_k = %5.5g'%sigma_PA_k)
    
    #cosmic variance only
    VAR_PA_k = (9 / 20 / np.pi) * sigma_PA_k**2 + sigma_P_k_21 * sigma_P_k_CO - sigma_PS_k**2
    VAR_PA_k_cv = (9 / 20 / np.pi) * sigma_PA_k**2 + sigma_P_k_21_cv * sigma_P_k_CO_cv - sigma_PS_k**2 #eqn 9 in Tan's draft
    sigma_PA_k, sigma_PA_k_cv = VAR_PA_k**0.5, VAR_PA_k_cv**0.5
    
    print('sigma_error = %5.5g'%sigma_PA_k, 'sigma_cv = %5.5g'%sigma_PA_k_cv)
    return sigma_PA_k, sigma_PA_k_cv

def significance_level(signal_noise_ratio_square_sum, free_degree):
    """
    Calculate the significance level for claiming a detection.
    Assume the square sum of signal-noise-ratio from multiple measurements follows a chi-square distribution,
    then get the p-value and convert it to the equivalent significance level of a normal distribution.
    -----------
    Parameters:
    -----------
    signal_noise_ratio_square_sum : the sum of the square of Signal to Noise ratio
    free_degree : the number of measurements minus 1
    --------
    Returns:
    --------
    significance : the significance level
    """
    p = chi2.cdf(signal_noise_ratio_square_sum, free_degree)
    significance = erfinv(p)
    if significance == float('inf'):
        p_comp = integrate.quad(lambda x: chi2.pdf(x, free_degree), signal_noise_ratio_square_sum, np.inf, epsrel = 1e-5)[0]
        significance = erfcinv(p_comp)
    return significance

def maximum_significance_level(kh_array, P_k, sigma_P_k):
    '''
    calculate the maximum significance level
    -----------
    Parameters:
    -----------
    kh_array : h/Mpc, the array of the index in k-space
    P_k : [muK]^2[Mpc/h]^3, the power spectrum
    sigma_P_k : [muK]^2[Mpc/h]^3, the error of the power spectrum
    --------
    significance : the maximum significance level
    kh_max : h/Mpc, the largest k within which the significance reaches the maximum
    '''
    free_degree = len(k_array) - 1
    S_over_N = [abs(Pk[i]) / sigma_P_k[i] for i in range(free_degree + 1)]
    square_sum = np.sum(np.array(S_over_N)**2)
    confidence = significance_level(square_sum, free_degree)
    kh_min = kh_array[0]
    kh_max = kh_array[free_degree]
    while free_degree > 1:
        square_sum -= S_over_N[free_degree]**2
        free_degree -= 1
        temp_confidence = significance_level(square_sum, free_degree)
        #print(square_sum, free_degree, temp_confidence)
        if (temp_confidence > confidence):
            confidence = temp_confidence
            kh_max = kh_array[free_degree]
    return confidence, kh_max
