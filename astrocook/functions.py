from .vars import *
from astropy import constants as ac
from scipy.special import wofz
#from lmfit.lineshapes import gaussian as gauss
import numpy as np

def _fadd(a, u):
    """ @brief Real part of the Faddeeva function Re(F)
    @param a First abstract variable
    @param u Second abstrac variable
    @return Re(F(a, u))
    """

    return np.real(wofz(u + 1j * a))


def adj_gauss(x, z, ampl, sigma, series='Ly_a'):
    model = np.ones(len(x))
    for t in series_d[series]:
        c = (1+z)*xem_d[t].to(au.nm).value
        model += ampl*np.exp(-(0.5 * (x-c) / sigma)**2)
    return model

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

def lines_voigt(x, z, N, b, btur, series='Ly_a'):
    """ @brief Voigt function (real part of the Faddeeva function, after a
    change of variables)

    @param x Wavelength domain (in nm)
    @param z Redshift
    @param N Column density (in cm^-2)
    @param b Doppler broadening (in km s^-1)
    @param btur Turbulent broadening (in km s^-1)
    @param series Series of ionic transition
    @param xem Wavelength of the line (in nm)
    @param tab Table with the Faddeeva function
    @return Voigt function over x
    """

    x = x * au.nm
    z = z * au.dimensionless_unscaled
    N = N / au.cm**2
    b = b * au.km/au.s
    btur = btur * au.km/au.s

    model = np.ones(len(x))
    for t in series_d[series]:
        xem = xem_d[t]
        xobs = xem*(1+z)
        fosc = fosc_d[t]
        gamma = gamma_d[t]/au.s
        b_qs = np.sqrt(b**2 + btur**2)
        atom = fosc * ac.e.esu**2 / (ac.m_e * ac.c)
        tau0 = np.sqrt(np.pi) * atom * N * xem / b_qs
        a = 0.25 * gamma * xem / (np.pi * b_qs)
        u = ac.c/b_qs * ((x/xobs).to(au.dimensionless_unscaled) - 1)

        model *= np.array(np.exp(-tau0.to(au.dimensionless_unscaled) \
                          * _fadd(a, u)))

    return model

def nfwhm_voigt(N=0, b=20, btur=0, series='Ly_a'):
    """ @brief FWHM of the Voigt function normalized to the speed of light.
    """

    nfwhm = []
    for t in series_d[series]:
        xem = xem_d[t]
        gamma = gamma_d[t]/au.s
        b_qs = np.sqrt(b**2 + btur**2)
        a = 0.25 * gamma * xem / (np.pi * b_qs)


        nfwhm.append(0.5346*a.value+np.sqrt(0.2166*a.value**2+1))
    ret = 2*b/ac.c.to(au.km/au.s).value
    print(nfwhm, ret)
    return ret



def psf_gauss(x, #center, resol):
              z, resol, series='Ly_a'):
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

    """
    sigma = center / resol * 4.246609001e-1
    psf = np.exp(-(0.5 * (x-center) / sigma)**2)
    psf[np.where(psf < 1e-4)] = 0.0
    ret = [np.array(psf)]#[c_min:c_max]]
    return ret
    """

    ret = []
    for t in series_d[series]:
        c = (1+z)*xem_d[t].value
        sigma = c / resol * 4.246609001e-1
        psf = np.exp(-(0.5 * (x-c) / sigma)**2)
        psf[np.where(psf < 1e-4)] = 0.0
        ret.append(psf)#[c_min:c_max]]
    return ret
