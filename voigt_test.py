from astrocook import *
from lmfit import Parameters
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import time

def main():



    # Input parameters
    #name = '../astrocook_data/B2126-15'
    name = 'B2126-15_part'
    #name = 'J0940_part'
    zem = 3.268

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Convolve the 1D spectrum with a gaussian filter
    sigma = 5.0
    conv = spec.convolve(gauss_sigma=sigma)

    # Find the lines in the 1D spectrum
    kappa = 5.0
    lines = conv.find_lines(mode='abs', kappa=kappa, diff='max')

    # Create a "voigt" object with the spectrum and the lines
    # Plot results
    fig = plt.figure(figsize = (10, 10))
    fig.canvas.set_window_title('Lines')
    fig.suptitle(name)

    #gs = gridspec.GridSpec(2, 3)
    gs = gridspec.GridSpec(2, 2) 
    ax_0 = fig.add_subplot(gs[0, :])
    ax_10 = fig.add_subplot(gs[1, :])
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

    for l in range(len(lines.x)):
    #for l in [0]:
        start_line = time.time()

        ax_0.scatter(lines.x[l], lines.y[l], s=100, color='g')

        voigt = Voigt(spec, lines)
        voigt._redchi = float('inf')
        #voigt.prep(l)
        #voigt.model_cont(l)
        #trasm_model = voigt.model_trasm(l)
        out = voigt.fit_auto(l)

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

        
        """
        ax_11.cla()
        ax_11.plot(voigt._x_ran, voigt._tau_rect, c='b')
        ax_11.plot(voigt._x_ran, voigt._tau_trasm, c='r', linestyle=':')
        ax_11.plot(voigt._x_ran, voigt._tau_norm0, c='y', linestyle=':')        
        ax_11.plot(voigt._x_ran, voigt._tau_fit, c='g')
        ax_11.plot(voigt._x_ran, voigt._tau_norm, c='y')        
        """
        
        plt.draw()
        plt.pause(0.001)

    print("Time to fit all the lines: %.0f seconds." \
          % (time.time() - start_all))
    plt.ioff()
    plt.show()
        
if __name__ == '__main__':
    main()
