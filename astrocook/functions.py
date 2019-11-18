from .vars import *
from astropy import constants as ac
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from scipy.special import wofz
#from lmfit.lineshapes import gaussian as gauss
from matplotlib import pyplot as plt
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

def convolve(data, psf):
    s = 0
    l = 0
    for i, k in enumerate(psf):
        s += l
        l = len(k)
        k_arr = k[np.where(k>0)]
        k_arr = k_arr/np.sum(k_arr)
        data_arr = data#[s:s+l]
        pad_l = len(k_arr)#*2
        pad = np.ones(pad_l)
        temp_arr = np.concatenate((pad*data_arr[0], data_arr, pad*data_arr[-1]))
        conv = np.convolve(temp_arr, k_arr, mode='valid')[pad_l//2+1:][:l]
        #plt.plot(range(len(conv)),temp_arr[1:-1], c='red', alpha=0.5)
        #plt.plot(range(len(k)),k)
        #plt.plot(range(len(conv)),conv, c='black', alpha=0.5)
        plt.show()
        if i == 0:
            ret = conv
        else:
            ret *= conv
    return ret

def detect_local_minima(arr):
    #https://stackoverflow.com/questions/3986345/how-to-find-the-local-minima-of-a-smooth-multidimensional-array-in-numpy-efficie
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    #
    # we create the mask of the background
    background = (arr==0)
    #
    # a little technicality: we must erode the background in order to
    # successfully subtract it from local_min, otherwise a line will
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    #
    # we obtain the final mask, containing only peaks,
    # by removing the background from the local_min mask
    detected_minima = local_min #- eroded_background
    #return np.where(detected_minima)
    return detected_minima

def lines_voigt(x, z, logN, b, btur, series='Ly_a'):
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

    #x = x[0] * au.nm
    x = x * au.nm
    z = z * au.dimensionless_unscaled
    N = 10**logN / au.cm**2
    b = b * au.km/au.s
    btur = btur * au.km/au.s
    model = np.ones(len(x))
    for t in series_d[series]:
        if series == 'unknown':
            xem = z*au.nm
            xobs = z*au.nm
        else:
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
        #model *= np.array(-tau0.to(au.dimensionless_unscaled) * _fadd(a, u)))

    return model

def parse(series):
    trans = []
    for s in series.split(','):
        if '_' in s:
            trans.append(s)
        else:
            for t in series_d[s]:
                trans.append(t)
    return trans

def psf_gauss(x, #center, resol):
              resol, reg=None):
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

    c = np.median(reg)
    sigma = c / resol * 4.246609001e-1
    psf = np.exp(-(0.5 * (x-c) / sigma)**2)
    #psf[np.where(psf < 1e-4)] = 0.0
    #psf = np.zeros(len(x))
    #psf[len(x)//2] = 1
    ret = [np.array(psf)]
    return ret

def running_mean(x, h=1):
    """ From https://stackoverflow.com/questions/13728392/moving-average-or-running-mean """

    n = 2*h+1
    cs = np.cumsum(np.insert(x, 0, 0))
    rm = (cs[n:] - cs[:-n]) / float(n)
    return np.concatenate((h*[rm[0]], rm, h*[rm[-1]]))
