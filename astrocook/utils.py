import numpy as np
from math import factorial
from astropy.constants import c
from scipy.special import wofz

def many_gauss(x, *p, mode='abs', cont = 1):
    """Sum of gaussian profiles
    
    Adapted from: 
    http://stackoverflow.com/questions/26902283/fit-multiple-gaussians-to-the-data-in-python
    """

    # Compute the profile
    sum = np.zeros_like(x)
    for i in range(0, len(p), 3):
        mean = p[i]
        std  = p[i + 1]
        ampl = p[i + 2]
        sum = sum + ampl * np.exp( -((x - mean) / std)**2)

    # Convert the profile in transmission, if is meant for absorption lines
    # (default)
    y = sum
    
    return y
    
def many_voigt(x, *p, constr_list=None, trans_list=None, transf='yes'):

    """
    Return the voigtoigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """

    trans_len = int(len(p) / 3)
    if trans_list is None:
        #trans_list = ['Lya'] * trans_len
        wave = [121.567] * trans_len
        oscs = [0.416] * trans_len
        damp = [6.265e8] * trans_len
    """
    wave = np.zeros(trans_len)
    oscs = np.zeros(trans_len)
    damp = np.zeros(trans_len)
    for j in range(0, trans_len):
        if trans_list[j] == 'Lya':
            wave[j] = 121.567
            oscs[j] = 0.416
            damp[j] = 6.265e8
        if trans_list[j] == 'CIV_1548':
            wave[j] = 154.8204
            oscs[j] = 0.1899
            damp[j] = 2.643e8
        if trans_list[j] == 'CIV_1550':
            wave[j] = 155.0781
            oscs[j] = 0.09475
            damp[j] = 2.628e8
    """
    
    #print(p)
    fact = 3600
    norm = [500.0 * fact, 14.0 * fact, 20.0 * fact] * int(len(p) / 3)
    #print(norm)
    if transf == 'yes':
        #print(p)
        p = norm * (np.arctan(p) + 0.0)
        #p = norm * np.array(p)
        #print(p)
            
    sum = np.zeros_like(x)
    for i in range(0, len(p), 3):
        j = int(i / 3)
        """
        w = p[i]
        log_N = p[i + 1] + 4.0
        b = p[i + 2] * 1e3
        u = c.value / (p[i + 2] * 1e3) * (1 - x / p[i])
        l_hwhm = 0.25 * damp[j] * wave[j] * 1.0e-9 / (np.pi * p[i + 2] * 1e3)
        g_hwhm = np.sqrt(2 * np.log(2))
        ampl = oscs[j] * np.power(10, p[i + 1] + 4.0) * 2.8179403267e-15 * np.sqrt(np.pi) \
               * c.value * wave[j] * 1.0e-9 / (p[i + 2] * 1e3)
        sum = sum + np.sqrt(np.pi) * np.real(wofz((u + 1j * l_hwhm) * g_hwhm)) \
              * g_hwhm / np.sqrt(np.pi) * ampl
        """

        u = c.value / (p[i + 2] * 1e3) * (1 - x / p[i])
        l_hwhm = 0.25 * damp[j] * wave[j] * 1.0e-12 / (np.pi * p[i + 2])
        g_hwhm = np.sqrt(2 * np.log(2))
        ampl = oscs[j] * np.power(10, p[i + 1]) * 2.8179403267 * np.sqrt(np.pi) \
               * c.value * wave[j] * 1.0e-23 / (p[i + 2])
        sum = sum + np.sqrt(np.pi) * np.real(wofz((u + 1j * l_hwhm) * g_hwhm)) \
              * g_hwhm / np.sqrt(np.pi) * ampl


    return np.exp(-sum)
    
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
       
    Adapted from:
    http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')