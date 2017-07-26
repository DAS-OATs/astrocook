from astrocook import *
import copy
from lmfit import Parameters
from matplotlib import gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import RectBivariateSpline, bisplrep, bisplev, interp2d
from scipy.special import wofz
import time
import sys

def main():



    # Input parameters
    #name = 'B2126-15_part2'
    #name = 'B2126-15_part3'    
    #name = 'J0940_part2'
    #name = 'J0940_CIV'
    name = 'J0940_CIV2'    
    zem = 3.268

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Convolve the 1D spectrum with a gaussian filter
    sigma = 4.0
    conv = spec.convolve(gauss_sigma=sigma)

    # Find the lines in the 1D spectrum
    kappa = 5.0
    lines = conv.find_lines(mode='abs', kappa=kappa, diff='max')

    # Find the continuum
    spec = Spec1DCont(spec)
    list = copy.deepcopy(lines)
    list.group()
    spec.rem_lines(list, size=5)
    spec.sg(spec.abs_rem, 101, 3)
    #print(spec._t)

    # Create a "voigt" object with the spectrum and the lines
    # Plot results
    spec_fig = plt.figure(figsize = (10, 2))
    spec_fig.canvas.set_window_title('Spectrum')
    spec_fig.suptitle(name)
    syst_fig = plt.figure(figsize = (6, 6))
    syst_fig.canvas.set_window_title('Lines')
    syst_fig.suptitle(name)

    spec_gs = gridspec.GridSpec(1, 1)
    spec_ax = spec_fig.add_subplot(spec_gs[:, :])
    spec_gs.tight_layout(spec_fig, rect=[0.01, 0.01, 1, 0.97])
    syst_gs = gridspec.GridSpec(1, 1)
    syst_ax = syst_fig.add_subplot(syst_gs[:, :])
    syst_gs.tight_layout(syst_fig, rect=[0.01, 0.01, 1, 0.97])

    spec_ax.plot(spec.x, spec.y, c='b')
    spec_ax.plot(conv.x, conv.y, c='r')
    spec_ax.plot(spec.x, spec.cont, c='r', linestyle=':')
    spec_ax.scatter(lines.x, lines.y)
    plt.ion()
   
    start_all = time.time()

    try:
        a_range = np.loadtxt('voigt_a_range.dat')
        u_range = np.loadtxt('voigt_u_range.dat')                
    except:
        all = Voigt(spec, lines)
        all.tabulate(pref='voigt_')
        a_range = np.loadtxt('voigt_a_range.dat')
        u_range = np.loadtxt('voigt_u_range.dat')                

    l = 0
    while l < len(lines.x):

        spec_ax.scatter(lines.x[l], lines.y[l], s=100, color='g')

        try:
            cont = voigt._y_cont
            redchi = voigt._out_redchi
            aic = voigt._out_aic            
        except:
            cont = None
            redchi = float('inf')
            aic = float('inf')
        """
        voigt = Voigt(spec, lines, chosen=l)
        print(voigt.t)
        """
        #"""
        try:
            lines_id = lines.t['SYST']
        except:
            lines_id = [['CIV_1548', 'CIV_1550'], ['CIV_1548', 'CIV_1550'], 
                        ['CIV_1550', 'CIV_1548'], ['CIV_1550', 'CIV_1548']]
        voigt = Voigt(spec, lines, chosen=l, syst=lines_id)
        #"""
        """
        lines_id = [['CIV_1548', 'CIV_1550'], ['CIV_1548', 'CIV_1550'], 
                    ['CIV_1548', 'CIV_1550'], ['CIV_1548', 'CIV_1550'], 
                    ['CIV_1548', 'CIV_1550'], ['CIV_1548', 'CIV_1550'], 
                    ['CIV_1548', 'CIV_1550'], ['CIV_1548', 'CIV_1550'], 
                    ['CIV_1550', 'CIV_1548'], ['CIV_1550', 'CIV_1548'], 
                    ['CIV_1550', 'CIV_1548'], ['CIV_1550', 'CIV_1548'], 
                    ['CIV_1550', 'CIV_1548'], ['CIV_1550', 'CIV_1548']]
        voigt = Voigt(spec, lines, chosen=l, syst=lines_id)
        """
        voigt.tab = np.loadtxt('voigt_tab.dat')
        voigt.splrep = interp2d(a_range, u_range, voigt.tab)

        try:
            if np.min(voigt.group(l)['XMIN']) \
               != np.min(voigt.group(l - 1)['XMIN']):
                cont = voigt.range(l)['CONT']
                redchi = float('inf')
                aic = float('inf')
        except:
            if cont is None:
                cont = voigt.range(l)['CONT']
                
        start_line = time.time()
        voigt.fit_auto(l, cont, redchi, aic, ax=syst_ax)
        
        spec_ax.plot(voigt._x_ran, voigt._y_cont0, c='y', linestyle=":")
        spec_ax.plot(voigt._x_ran, voigt._y_cont, c='y')
        try:
            spec_ax.plot(voigt._x_ran, voigt._y_fit, c='g')
        except:
            pass

        syst_ax.cla()
        comp = voigt.group(l)
        for x in comp['X']:
            syst_ax.axvline(x=x, ymin=0.75, ymax=0.95, color='lightgray')
        syst_ax.plot(voigt._x_ran, voigt._y_ran / voigt._y_cont, c='b')
        syst_ax.plot(voigt._x_ran, voigt._dy_ran / voigt._y_cont, c='b',
                   linestyle=':')        
        syst_ax.plot(voigt._x_ran, -voigt._dy_ran / voigt._y_cont, c='b',
                   linestyle=':')
        syst_ax.plot(voigt._x_ran, voigt._y_trasm / voigt._y_cont, c='r',
                   linestyle=':')
        syst_ax.plot(voigt._x_ran, voigt._y_cont0 / voigt._y_cont, c='y',
                   linestyle=':')
        try:
            syst_ax.plot(voigt._x_ran, voigt._y_fit / voigt._y_cont, c='g')
            syst_ax.plot(voigt._x_ran, voigt._y_cont / voigt._y_cont, c='y')
            syst_ax.plot(voigt._x_ran, voigt._y_resid / voigt._y_cont, c='g',
                       linestyle='--')
        except:
            pass
        plt.draw()
        plt.pause(0.1)
        
        lines._t = copy.deepcopy(voigt._t)
        print(lines._t)
        l += 1       
        print("Time to process line %2i of (at least) %3i: %.0f seconds." \
              % (l + 1, len(lines.x), time.time() - start_line))

        
    print("Time to process the spectrum: %.0f seconds." \
          % (time.time() - start_all))
    plt.ioff()
    plt.show()
        
if __name__ == '__main__':
    main()
