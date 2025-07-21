import numpy as np
import matplotlib.pyplot as plt

c = 299792458.0 # m/s

def n_YAG(wl):
    """
    Returns the Sellmeier equation for YAG
    Args:
        wl: The wavelength to compute n at (um)
    returns:
        n: n computes at wl based on the Sellmeier equation (unitless)
    """
    # From Zelmon et al. see refractiveindex.com
    B1, C1 = 2.282, 0.01185
    B2, C2 = 3.27644, 282.734

    wl2 = wl**2
    n = np.sqrt(1
                + B1*wl2/(wl2 - C1)
                + B2*wl2/(wl2 - C2)
                )
    return n

def n_Si(wl):
    """
    Returns the Sellmeier equation for Si
    Args:
        wl: The wavelength to compute n at (um)
    returns:
        n: n computes at wl based on the Sellmeier equation (unitless)
    """
    # From Salzberg and Villa see refractiveindex.com
    B1, C1 = 10.6684293, 0.301516485
    B2, C2 = 0.0030434748, 1.13475115
    B3, C3 = 1.54133408, 1104

    wl2 = wl**2
    n = np.sqrt(1
                + B1*wl2/(wl2 - C1**2)
                + B2*wl2/(wl2 - C2**2)
                + B3*wl2/(wl2 - C3**2)
                )
    return n

def n_CaF2(wl):
    """
    Returns the Sellmeier equation for Si
    Args:
        wl: The wavelength to compute n at (um)
    returns:
        n: n computes at wl based on the Sellmeier equation (unitless)
    """
    # From Lisee refractiveindex.com
    B1, C1 = 0.69913, 0.09374
    B2, C2 = 0.11994, 21.18
    B3, C3 = 4.35181, 38.46

    wl2 = wl**2
    n = np.sqrt(1 + 0.33973
                + B1*wl2/(wl2 - C1**2)
                + B2*wl2/(wl2 - C2**2)
                + B3*wl2/(wl2 - C3**2)
                )
    return n

def wl_to_omega(wl):
    """Convert wavelength (m) to angular frequency (rad/s)"""
    return 2 * np.pi * c / wl

def omega_to_wl(omega):
    """Conver angular frequency (rad/s) to wavelength (m)"""
    return 2 * np.pi * c / omega

def omega_grid(N, center_wl, fractional_range):
    """
    Returns an angular frequency grid centered at the center wl. fractional range
    determines what fraction of the center wavelength +/- the range is.

    Args:
        N: number of grid points must be odd
        center_wl: center wavelength of grid
        fractional_range: fraction of the center wavelength +/-, determines grid bounds
    returns:
        omega_grid: omega grid centered at center_wl   
    """

    if (N % 2) == 0:
        raise ValueError("N must be odd")
    if fractional_range >= 1:
        raise ValueError("Fractional range must be less than 1")

    omega0 = wl_to_omega(center_wl)
    omega_range = omega0 * fractional_range
    omega_grid = omega0 + np.linspace(-omega_range, omega_range, N)
    return omega_grid

def propagation_const(omega, sellmeier_equation):
    """
    Computes the propagation constant beta(omega) for some sellmeier_equation at some wavelength

    Args:
        omega: omega axis centered around the frequency of interest
        sellmeier_equation: sellmeier equation centered at the frequency (wavelength) of interest
    returns:
        beta: returns the propagation cosntant as a function of angular freq
    """
    if type(omega) != np.ndarray and type(omega) != np.array:
        raise TypeError("Omega axis must be of type np.ndarray or np.array")
    if type(sellmeier_equation) != np.ndarray and type(sellmeier_equation) != np.array:
        raise TypeError("Sellmeier equation must be of type np.ndarray or np.array")
    return (sellmeier_equation * omega) / c

def dispersion_coeff(beta0, omega):
    """
    Calculates the second, third, and fourth dispersion coefficients. keep in mind your beta 0
    should have been made with the same omega axis that you pass to this function.

    Args:
        beta0: an array containing propagation constant values (1/m)
        omega: angular frequency axis centered around the frequency of interest (rad/s)
    returns:
        beta2, beta3, beta4: returns dispersion coefficients of order 2,3,4
        (s^2/m)  (s^3/m)  (s^4/m)
    
    """
    if len(omega) % 2 == 0:
        raise ValueError("Length of omega axis must be odd")
    if len(beta0) % 2 == 0:
        raise ValueError("Length of beta0 axis must be odd")
    if len(beta0) != len(omega):
        raise ValueError("beta0 and omega axis must be the same length")
    
    d1 = np.gradient(beta0, omega)
    d2 = np.gradient(d1, omega)
    d3 = np.gradient(d2, omega)
    d4 = np.gradient(d3, omega)

    N = len(omega)

    beta2 = d2[N//2]
    beta3 = d3[N//2]
    beta4 = d4[N//2]

    return beta2, beta3, beta4

def fwhm(x, y, height = 0.5):
    """
    Gives the full-width at half-maximum of data in a numpy array pair.

    :param x: The x-values (e.g. the scale of the vector; has the same units as the return value)
    :type x: np.ndarray

    :param y: The y-values (to which half-maximum refers)
    :type y: np.ndarray

    :param height: Instead of half-max, can optionally return height*max (e.g. default 0.5)
    :type height: float

    :return: The full-width at half-max (units of x)
    :rtype: float
    """
    heightLevel = np.max(y) * height
    indexMax = np.argmax(y)
    y = np.roll(y, - indexMax + int(np.shape(y)[0]/2),axis=0)
    indexMax = np.argmax(y)
    xLower = np.interp(heightLevel, y[:indexMax], x[:indexMax])
    xUpper = np.interp(heightLevel, np.flip(y[indexMax:]), np.flip(x[indexMax:]))
    return xUpper - xLower

def map_beam_radius(beam_waist, rayleigh_length, max_z, precision, verbose = True):
    """
    Calculates the beam waist of the described guassian beam for some extent of z. returns a list of these
    values and has the option to plot the function.
    
    Args:
        beam_waist: guassian beam waist, most narrow radius of beam (m)
        rayleigh_length: desired rayleigh length of beam (m)
        max_z: extent to which one wants their array created to (m)
        precision: what spatial increment to compute the beam radius at, starts at 0 (m)
        verbose: whether to graph the results or not

    returns:
        tuple: A tuple containing the z_axis used for computation as well as the beam radius values
    
    """
    if (max_z < precision):
        raise ValueError("Maximum Z value must be greater than the desired precision")
    
    z_axis = np.arange(0.0, max_z, precision)
    beam_radius = np.sqrt( beam_waist**2 * (1 + (z_axis / rayleigh_length)**2) )

    if verbose:
        plt.figure(figsize=(10,5))
        plt.plot(z_axis, beam_radius, label=r'Beam Radius')
        plt.xlabel("Z Position (m)")
        plt.ylabel("Beam Radius (m)")
        plt.title(" Beam Radius vs Z ")
        plt.grid()
        if max_z > rayleigh_length:
            plt.axvline(x=rayleigh_length, color='red', linestyle='--', label='$Z_{R}$')
        plt.legend()
        plt.show()

    return (z_axis, beam_radius)

def dB(num):
    return 10*np.log10(np.abs(num)**2)

def freq_to_wl(freq_THz):
    """Convert frequency in THz to wavelength in um"""
    c = 299792458  # Speed of light in m/s
    return c / (freq_THz * 1e12) * 1e6  # Convert to m

def eff_nonlin(n2, wl, radius):
    """
    Returns the effective nonlinearity in (1/Wm) given the central wl, nonlinear refractive index,
    and fiber/beam radius.

    Args:
        n2: nonlinear refractive index (m^2/W)
        wl: central wavelength (m)
        radius: effective radius (m)
    return:
        gamma: Effective nonlinearity (1/Wm)
    """
    gamma = (2 * np.pi * n2) / (wl * np.pi * radius**2)
    return gamma

def calc_radius(pulse_energy, pulse_duration, peak_intensity):
    """
    Returns the beam radius when given the pule energy, duration, and peak intensity

    Args:
        pulse_energy: pulse energy (J)
        pulse_duration: pulse duration (s)
        peak_intesnity: peak pulse intensity (W/m^2)
    Returns:
        radius: calculated beam radius (m)   
    """
    radius = np.sqrt(
        (1.88 * pulse_energy) / (np.pi * peak_intensity * pulse_duration)
    )
    return radius


