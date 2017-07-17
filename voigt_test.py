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

def main():



    # Input parameters
    #name = '../astrocook_data/B2126-15'
    name = 'B2126-15_part2'
    #name = 'B2126-15_part3'    
    #name = 'J0940_part2'
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
    #print(cont._t)

    # Create a "voigt" object with the spectrum and the lines
    # Plot results
    fig = plt.figure(figsize = (10, 7))
    fig.canvas.set_window_title('Lines')
    fig.suptitle(name)

    gs = gridspec.GridSpec(2, 3)
    ax_0 = fig.add_subplot(gs[0, :])
    ax_10 = fig.add_subplot(gs[1, :])
    gs.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])

    ax_0.plot(spec.x, spec.y, c='b')
    ax_0.plot(conv.x, conv.y, c='r')
    ax_0.plot(spec.x, spec.abs_rem, c='r', linestyle=':')
    ax_0.scatter(lines.x, lines.y)
    plt.ion()
    #plt.show()
   
    start_all = time.time()
    #for l in range(len(lines.x)):

    try:
        a_range = np.loadtxt('voigt_a_range.dat')
        u_range = np.loadtxt('voigt_u_range.dat')                
    except:
        all = Voigt(spec, lines)
        all.tabulate(pref='voigt_')
        a_range = np.loadtxt('voigt_a_range.dat')
        u_range = np.loadtxt('voigt_u_range.dat')                

    #for l in range(len(lines.x)):
    #for l in [0]:
    #for l in [len(lines.x) - 1]:
    l = 0
    while l < len(lines.x):

    
        ax_0.scatter(lines.x[l], lines.y[l], s=100, color='g')

        voigt = Voigt(spec, lines, chosen=l)        
        voigt.tab = np.loadtxt('voigt_tab.dat')
        voigt.splrep = interp2d(a_range, u_range, voigt.tab)
        
        try:
            if np.min(voigt.group(l)['XMIN']) \
               == np.min(voigt.group(l - 1)['XMIN']):
                cont = voigt._y_cont
                redchi = voigt._out_redchi
                aic = voigt._out_aic
            else:
                cont = voigt.range(l)['ABS_REM']
                redchi = float('inf')
                aic = float('inf')
        except:
            cont = voigt.range(l)['ABS_REM']
            redchi = float('inf')
            aic = float('inf')

        
        start_line = time.time()
        #print(cont)
        voigt.fit_auto(l, cont, redchi, aic, ax=ax_10)
        #print(voigt.group(l))
        print("Time to process line %2i of %3i (at least): %.0f seconds." \
              % (l + 1, len(lines.x), time.time() - start_line))

        ax_0.plot(voigt._x_ran, voigt._y_cont0, c='y', linestyle=":")
        ax_0.plot(voigt._x_ran, voigt._y_cont, c='y')
        try:
            ax_0.plot(voigt._x_ran, voigt._y_fit * voigt._y_slope, c='g')
        except:
            pass

        ax_10.cla()
        comp = voigt.group(l)
        for x in comp['X']:
            ax_10.axvline(x=x, ymin=0.75, ymax=0.95, color='lightgray')
        ax_10.plot(voigt._x_ran, voigt._y_rect, c='b')
        ax_10.plot(voigt._x_ran, voigt._dy_rect, c='b', linestyle=':')        
        ax_10.plot(voigt._x_ran, -voigt._dy_rect, c='b', linestyle=':')
        ax_10.plot(voigt._x_ran, voigt._y_trasm, c='r', linestyle=':')
        ax_10.plot(voigt._x_ran, voigt._y_norm0, c='y', linestyle=':')
        try:
            ax_10.plot(voigt._x_ran, voigt._y_fit, c='g')
            ax_10.plot(voigt._x_ran, voigt._y_norm, c='y')
            ax_10.plot(voigt._x_ran, voigt._y_resid, c='g', linestyle='--')
        except:
            pass
        plt.draw()
        plt.pause(0.1)
        
        lines._t = copy.deepcopy(voigt._t)
        l += 1       

    print("Time to process the spectrum: %.0f seconds." \
          % (time.time() - start_all))
    plt.ioff()
    plt.show()
        
if __name__ == '__main__':
    main()
