from astrocook import List, ListReader, ListSyst, Spec1DCont, Spec1DReader
from astrocook.utils import many_voigt, many_gauss
from astropy import units as u
from astropy.constants import c
from astropy.table import vstack
import copy
from cycler import cycler
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import sys

import warnings
warnings.filterwarnings("ignore")

def in1d_tol(a, b, tol):
    """From http://stackoverflow.com/questions/23520952/python-change-format-of-np-array-or-allow-tolerance-in-in1d-function
    """
    
    d = np.abs(a - b[:, np.newaxis])
    return np.any(d<=tol, axis=0)


def main():
    """Find the absorption systems in a normalized QSO spectrum"""
    
    name = str(sys.argv[1])
    zem = str(sys.argv[2])
    
    fig = plt.figure(figsize = (12, 6))
    fig.canvas.set_window_title('Absorption systems')
    fig.suptitle(name)
    gs = gridspec.GridSpec(3, 1) 
    ax_0 = fig.add_subplot(gs[0])
    ax_1 = fig.add_subplot(gs[1])
    ax_2 = fig.add_subplot(gs[2])
    ax_0.set_prop_cycle(cycler('color', ['blue', 'red', 'green']))
    ax_1.set_prop_cycle(cycler('color', ['blue', 'cyan', 'lightgreen']))
    ax_2.set_prop_cycle(cycler('color', ['white', 'cyan', 'lightgreen']))
    #ax_0.set_ylim([-0.5, 1.5])   
    gs.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])



#    try:

    spec_read = Spec1DReader()
    spec = spec_read.cont(name + '_spec_cont.fits')
    list_read = ListReader()
    list_abs = list_read.gen(name + '_list_abs.fits')
    norm = copy.deepcopy(spec)
    norm.y = norm.y / norm.cont
    norm.dy = norm.dy / norm.cont
    print("Normalized spectrum loaded.")

    print("Shifting spectrum...")
    norm.todo_convert_logx(xunit=u.km/u.s)
    shift_l = copy.deepcopy(norm)
    shift_r = copy.deepcopy(norm)
    trans_name = 'CIV'
    trans_wave = [154.8204, 155.0781] * u.nm
    shift_l.x = norm.x - c * (1 - trans_wave[0] / trans_wave[1])
    shift_r.x = norm.x + c * (1 - trans_wave[0] / trans_wave[1])
    shift_l.y = np.interp(norm.x, shift_l.x, shift_l.y)
    shift_r.y = np.interp(norm.x, shift_r.x, shift_r.y)
    shift_l.x = norm.x
    shift_r.x = norm.x
    norm.todo_convert_logx(xunit=u.nm)

    print("Smoothing shifted spectra...")
    # Higher values of gauss_sigma decrease the number of false positives
    conv = norm.convolve(gauss_sigma=20.0)
    conv_l = shift_l.convolve(gauss_sigma=20.0)
    conv_r = shift_r.convolve(gauss_sigma=20.0)

    print("Correlating shifted spectra...")
    corr_l = copy.deepcopy(conv)
    corr_r = copy.deepcopy(conv)
    corr_l.y = 1 - np.sqrt((1 - conv.y.value) * (1 - conv_l.y.value))
    corr_r.y = 1 - np.sqrt((1 - conv.y.value) * (1 - conv_r.y.value))

    conv_l.todo_convert_logx(xunit=u.nm)
    conv_r.todo_convert_logx(xunit=u.nm)

    print("Finding " + trans_name + " candidates...")

    # Find candidates
    #conv_l = corr_l.convolve(gauss_sigma=5.0)
    corr_l._t = corr_l._t[corr_l.x.value > (1.0 + float(zem)) * 121.567]
    corr_r._t = corr_r._t[corr_r.x.value > (1.0 + float(zem)) * 121.567]

    # Higher valeus of kappa decrease the number of false positives
    cand_l = corr_l.find_lines(kappa=10.0)
    cand_r = corr_r.find_lines(kappa=10.0)
    z_l = cand_l.x.value / trans_wave[0].value
    z_r = cand_r.x.value / trans_wave[1].value
    
    # Select outliers with a doublet companion
    z_l_coin = in1d_tol(z_l, z_r, 1e-4)
    z_r_coin = in1d_tol(z_r, z_l, 1e-4)
    cand_l._t = cand_l._t[z_l_coin]
    cand_r._t = cand_r._t[z_r_coin]
    print(" %i candidates found." % len(z_l_coin))
#    cand = List()
#    cand._t = vstack([cand_l.t, cand_r.t])
    
    # Create a line list of absorption systems, backtracing candidates in the
    # normalized spectrum to define Y values
    syst_l = ListSyst(cand_l, syst=range(1, len(cand_l.t) + 1), 
                      name=['CIV_1548'] * len(cand_l.t), z=z_l[z_l_coin])
    norm_l_coin = np.in1d(norm.x.value, syst_l.x.value)
    syst_l.y = norm.y[norm_l_coin]
    syst_r = ListSyst(cand_r, syst=range(1, len(cand_r.t) + 1), 
                      name=['CIV_1550'] * len(cand_r.t), z=z_r[z_r_coin]) 
    norm_r_coin = np.in1d(norm.x.value, syst_r.x.value)
    syst_r.y = norm.y[norm_r_coin]

    syst = copy.deepcopy(syst_l)
    syst._t = vstack([syst_l.t, syst_r.t])

    syst.group()
    print(syst.t)
    """
    spec_cand = in1d_tol(list_abs.x.value, pre_x.value, 0.01)
    spec_cand_x = list_abs.x[spec_cand]
    spec_cand_y = list_abs.y[spec_cand]
    """

    # Plot results
    ax_0.plot(norm.x, norm.y)
    ax_0.scatter(syst.x, syst.y, s=200, c='black', marker='+')
    ax_0.set_ylabel("Y / CONT")
    ax_1.plot(conv.x, conv.y)
    ax_1.plot(conv_l.x, conv_l.y)
    ax_1.plot(conv_r.x, conv_r.y)
    ax_1.set_ylabel("Y / CONT")
    ax_2.plot(conv_l.x, conv_l.y)
    ax_2.plot(corr_l.x, corr_l.y)
    ax_2.plot(corr_r.x, corr_r.y)
    ax_2.scatter(cand_l.x, cand_l.y, s=200, c='blue', marker='+')
    ax_2.scatter(cand_r.x, cand_r.y, s=200, c='green', marker='+')
    #ax_2.scatter(list_l.x, list_l.y, s=200, c='black', marker='+')
    #ax_2.scatter(list_r.x, list_r.y, s=200, c='black', marker='+')
    ax_2.set_ylabel("Y / CONT")
    plt.show()

    print("Success!")


    #except:
    #    print("Normalized spectrum not found! Try running qso_cont.py first.")

if __name__ == '__main__':
    main()