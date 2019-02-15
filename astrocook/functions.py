from .vars import *
from scipy.special import wofz
import numpy as np

def fadd(a, u):
    """ @brief Real part of the Faddeeva function Re(F)
    @param a First abstract variable
    @param u Second abstrac variable
    @return Re(F(a, u))
    """

    return np.real(wofz(u + 1j * a))

def convolve(data, kernel):
    s = 0
    l = 0
    ret = np.array([])
    for k in kernel:
        s += l
        l = len(k)
        k_arr = k[np.where(k>0)]
        k_arr = k_arr/np.sum(k_arr)
        data_arr = data[s:s+l]
        pad_l = len(k_arr)#*2
        pad = np.ones(pad_l)
        temp_arr = np.concatenate((pad*data_arr[0], data_arr, pad*data_arr[-1]))
        conv = np.convolve(temp_arr, k_arr, mode='valid')[pad_l//2+1:][:l]
        #conv = np.convolve(data_arr, k_arr, mode='valid')[pad_l//2+1:][:l]
        ret = np.append(ret, conv)
    return ret


def gaussian(x,
             #c_min, c_max,
             center, resol):
    """ @brief Gaussian PSF

    The function returns a gaussian array for each element of a selected region
    in the wavelength domain

    @param x Wavelength domain (in nm)
    @param c_min Starting pixel of the region
    @param c_max Ending pixel of the region
    @param center Center wavelength of the region
    @param resol Resolution
    @return Gaussian PSF over x
    """

    sigma = center / resol * 4.246609001e-1
    psf = np.exp(-(0.5 * (x-center) / sigma)**2)
    psf[np.where(psf < 1e-4)] = 0.0
    ret = [np.array(psf)]#[c_min:c_max]]
    return ret


def voigt(x, z, N, b, btur, trans='Ly_a', xem=0.0, tab=None):
    """ @brief Voigt function (real part of the Faddeeva function, after a
    change of variables)

    @param x Wavelength domain (in nm)
    @param z Redshift
    @param N Column density (in cm^-2)
    @param b Doppler broadening (in km s^-1)
    @param btur Turbulent broadening (in km s^-1)
    @param trans Ionic transition
    @param xem Wavelength of the line (in nm)
    @param tab Table with the Faddeeva function
    @return Voigt function over x
    """

    x = x * au.nm
    z = z * au.dimensionless_unscaled
    N = N / au.cm**2
    b = b * au.km/au.s
    btur = btur * au.km/au.s


    if xem < 1e8:
        xem = xem_d[trans]#.to(au.m)
        xobs = xem*(1+z)
    else:  # For unknown species, Ly-alpha rest-frame wavelength is assumed
        print("xem", xem)
        xobs = xem
        xem = xem_d['Ly_a']
    fosc = fosc_d[trans]
    gamma = gamma_d[trans]/au.s
    b_qs = np.sqrt(b**2 + btur**2)
    atom = fosc * e.esu**2 / (m_e * c)
    tau0 = np.sqrt(np.pi) * atom * N * xem / b_qs
    a = 0.25 * gamma * xem / (np.pi * b_qs)
    u = c / b_qs * ((x / xobs).to(au.dimensionless_unscaled) - 1)

    model = np.exp(-tau0.to(au.dimensionless_unscaled) * fadd(a, u)) #* yunit
    return np.array(model)
