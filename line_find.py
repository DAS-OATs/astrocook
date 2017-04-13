from astrocook.spec_1d import Spec1D
from astrocook import Spec1DReader
from astropy import units as u
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import argrelmin, argrelmax
import sys

def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

def roll_wind(a, window):
    """ Define a rolling window 
    
    Source: http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html"""

    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def main():
    """Detect and fit absorption lines
    
    Algorithm:
        1. Undersample the spectrum to increase SNR;
        2. Find the extrema of the undersampled spectrum;
        3. Estimate the equivalent width of the absorption features;
        4. Select as lines the features with equivalent width above threshold.

    The equivalent width of absorption features is estimated by scanning the 
    list of extrema and computing the area of the triangles obtained by 
    connecting each local minimum to its closest significant maxima. Significant 
    maxima are selected among local maxima with a moving window of given width.
    """

    r = Spec1DReader()
    s = r.uves(str(sys.argv[1]))

    print("Detecting absorption lines...")
    
    # Undersample the spectrum to increase SNR
    try:
        bin = float(sys.argv[2])
    except:
        bin = 0.005
    x_rebin = np.trunc(s.x / bin)
    grouped = s.t.group_by(x_rebin)
    binned = grouped.groups.aggregate(np.mean)
    #s_b = r.table(binned)  
    s_b = s.from_table(binned) 
                 
    # Find the extrema of the undersampled spectrum           
    min = np.hstack(argrelmin(s_b.y))
    max = np.hstack(argrelmax(s_b.y))
    extr = np.sort(np.append(min, max))
    
    # Estimate the equivalent width of the absorption features
    try:
        width = int(sys.argv[3])
    except:
        width = 5
    x_win = roll_wind(s_b.x[extr], width)
    y_win = roll_wind(s_b.y[extr], width) 
    x_abs = x_win[np.arange(len(x_win)), np.argmax(y_win, 1)]
    y_abs = np.amax(y_win, 1)
    x_abs_l = s_b.x[extr[width - 1 : -width + 1]].value - x_abs[:-width + 1]
    x_abs_r = x_abs[width - 1:] - s_b.x[extr[width - 1 : -width + 1]].value
    y_abs_l = y_abs[:-width + 1] - s_b.y[extr[width - 1 : -width + 1]].value
    y_abs_r = y_abs[width - 1:] - s_b.y[extr[width - 1 : -width + 1]].value
    ew_abs = (x_abs_l * y_abs_r + x_abs_r * y_abs_l) \
        / (y_abs[:-width + 1] + y_abs[width - 1:])
    
    # Select as lines the features with equivalent width above threshold 
    try:
        kappa = float(sys.argv[4])
    except:
        kappa = 5.0
    clean_abs = sigma_clip(ew_abs)
    sel_idx = np.where(ew_abs > kappa * np.std(clean_abs))[0]
    ew_sel = ew_abs[sel_idx]
    x_sel_l = x_abs_l[sel_idx]
    x_sel_r = x_abs_r[sel_idx]
    extr_sel = extr[sel_idx + width - 1]
    extr_sort = np.argsort(extr_sel)
  
    # Create the list of lines (TBC)
    line = np.setdiff1d(extr_sel, max)
    line_pos = np.searchsorted(extr_sel[extr_sort], line)
    line_idx = extr_sort[line_pos]
    line_x = s_b.x[line]
    line_xmin = line_x - x_sel_l[line_idx] * line_x.unit
    line_xmax = line_x + x_sel_r[line_idx] * line_x.unit
    line_y = s_b.y[line]
    line_ew = ew_sel[line_idx] * line_x.unit

    print(" " + str(len(line_x)) + " lines detected.")
    print("Fitting absorption lines...")

    # Fit a gaussian profile
    fit_x = np.array([])
    fit_y = np.array([])
    fit_ew = np.array([])
    try:
        n_line = int(sys.argv[5])
    except:
        n_line = len(line_x)
    for l in range(0, len(line_x[0:n_line])):
        region = s.t[np.logical_and(s.x > line_xmin[l], s.x < line_xmax[l])]
        #s_r = r.table(region)
        s_r = s.from_table(region)
        print(s_r.t)
        print(len(s_r.t))
        prof_x = s_r.x.value
        prof_y = np.amax(s_r.y.value) - s_r.y.value
        ew_prof = line_ew[l].value
        fwhm_prof = ew_prof * np.amax(s_r.y.value) / np.amax(prof_y)
        p0 = [np.amax(s_r.y.value) - line_y[l].value, line_x[l].value, fwhm_prof / 2.355]
        print(p0)
        coeff, var_matrix = curve_fit(gauss, prof_x, prof_y, p0 = p0)
        fit_y_l = np.amax(s_r.y.value) - gauss(prof_x, *coeff)
        fit_x = np.append(fit_x, prof_x)
        fit_y = np.append(fit_y, fit_y_l)
        fit_ew = np.append(fit_ew, coeff[2] * 2.355 * np.amax(prof_y) \
            / np.amax(s_r.y.value))

    fit_sort = np.argsort(fit_x)
    fit_x = fit_x[fit_sort]
    fit_y = fit_y[fit_sort]

    print(" " + str(len(line_x[0:n_line])) + " lines fitted.")

    # Plot results
    fig = plt.figure(figsize = (15, 5))
    plt.plot(s.x, s.y, 'b-', s_b.x, s_b.y, 'r-')
    plt.plot(s_b.x[extr[width - 1 : -width + 1]], ew_abs * 1e3, 'g-')
    plt.plot(fit_x, fit_y, 'g-')
    plt.scatter(line_x, line_y, facecolors='none', edgecolors='r')
    plt.scatter(line_x, line_ew * 1e3, facecolors='none', edgecolors='g')
    plt.scatter(line_x[0:n_line], fit_ew * 1e3, facecolors='r', edgecolors='r')
    plt.show()
    """
    fig = plt.figure(figsize = (5, 5))
    plt.scatter(line_ew, fit_ew, color = 'r')
    plt.show()
    """  
if __name__ == '__main__':
    main()