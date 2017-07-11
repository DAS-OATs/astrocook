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
    #name = 'B2126-15_part2'
    name = 'B2126-15_part3'    
    #name = 'J0940_part2'
    zem = 3.268

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Convolve the 1D spectrum with a gaussian filter
    sigma = 10.0
    conv = spec.convolve(gauss_sigma=sigma)

    # Find the lines in the 1D spectrum
    kappa = 5.0
    lines = conv.find_lines(mode='abs', kappa=kappa, diff='max')

    # Create a "voigt" object with the spectrum and the lines
    # Plot results
    fig = plt.figure(figsize = (10, 10))
    fig.canvas.set_window_title('Lines')
    fig.suptitle(name)

    gs = gridspec.GridSpec(2, 3)
    ax_0 = fig.add_subplot(gs[0, :])
    ax_10 = fig.add_subplot(gs[1, :])
    #ax_10 = fig.add_subplot(gs[1, 0])
    #ax_11 = fig.add_subplot(gs[1, 1], projection='3d')
    #ax_11 = fig.add_subplot(gs[1, 1])    
    #ax_12 = fig.add_subplot(gs[1, 2])
    gs.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])

    ax_0.plot(spec.x, spec.y, c='b')
    ax_0.plot(conv.x, conv.y, c='r')
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
    for l in [len(lines.x) - 1]:


        ax_0.scatter(lines.x[l], lines.y[l], s=100, color='g')

        voigt = Voigt(spec, lines, chosen=l)
        voigt._redchi = float('inf')
        #voigt.prep(l)
        #voigt.model_cont(l)
        #trasm_model = voigt.model_trasm(l)
        #voigt.tab = copy.deepcopy(all.tab)
        #voigt.splrep = copy.deepcopy(all.splrep)
        voigt.tab = np.loadtxt('voigt_tab.dat')
        voigt.splrep = interp2d(a_range, u_range, voigt.tab)
        start_line = time.time()
        out = voigt.fit_auto(l, ax=ax_10)

        print("Time to fit the line: %.0f seconds." \
              % (time.time() - start_line))

        ax_0.plot(voigt._x_ran, voigt._y_cont0, c='y', linestyle=":")
        ax_0.plot(voigt._x_ran, voigt._y_cont, c='y')
        ax_0.plot(voigt._x_ran, voigt._y_fit * voigt._y_slope, c='g')

        ax_10.cla()
        comp = voigt.group(l)
        for x in comp['X']:
            ax_10.axvline(x=x, ymin=0.75, ymax=0.95, color='lightgray')
        ax_10.plot(voigt._x_ran, voigt._y_rect, c='b')
        ax_10.plot(voigt._x_ran, voigt._dy_rect, c='b', linestyle=':')        
        ax_10.plot(voigt._x_ran, -voigt._dy_rect, c='b', linestyle=':')
        ax_10.plot(voigt._x_ran, voigt._y_trasm, c='r', linestyle=':')
        ax_10.plot(voigt._x_ran, voigt._y_norm0, c='y', linestyle=':')
        ax_10.plot(voigt._x_ran, voigt._y_fit, c='g')
        ax_10.plot(voigt._x_ran, voigt._y_norm, c='y')
        ax_10.plot(voigt._x_ran, voigt._y_resid, c='g', linestyle='--')
        #ax_10.scatter(voigt.group(l)['X'], voigt.group(l)['Y'], s=100,
        #              color='g') 
        #ax_10.scatter(lines.x[l], lines.y[l], s=100,
        #              color='g') 

        #X, Y = np.meshgrid(voigt.pos_range, voigt.sigma_range)
        #Z = voigt.splrep(voigt.pos_range, voigt.sigma_range)
        #Z1 = voigt.tab
        #ax_11.plot_surface(X, Y, Z)
        #ax_11.plot_surface(X, Y, Z1)        
        """
        ax_11.cla()
        ax_11.plot(voigt._x_ran, voigt._tau_rect, c='b')
        ax_11.plot(voigt._x_ran, voigt._tau_trasm, c='r', linestyle=':')
        ax_11.plot(voigt._x_ran, voigt._tau_norm0, c='y', linestyle=':')        
        ax_11.plot(voigt._x_ran, voigt._tau_fit, c='g')
        ax_11.plot(voigt._x_ran, voigt._tau_norm, c='y')        
        """
        

    print("Time to fit all the lines: %.0f seconds." \
          % (time.time() - start_all))
    plt.ioff()
    plt.show()
        
if __name__ == '__main__':
    main()
