from astrocook import *
from matplotlib import gridspec
import matplotlib.pyplot as plt
import time

def main():

    start = time.time()

    # Input parameters
    #name = '../astrocook_data/B2126-15'
    name = 'B2126-15_part'
    zem = 3.268

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Convolve the 1D spectrum with a gaussian filter
    sigma = 2.0
    conv = spec.convolve(gauss_sigma=sigma)

    # Find the lines in the 1D spectrum
    kappa = 4.0
    lines = conv.find_lines(mode='abs', kappa=kappa, diff='max')

    # Create a "fit" object with the spectrum and the lines
    fit = Fit(spec, lines)
    print(fit.t)
    #print(fit.group())
    #print(fit.comp())
    #print(fit.range())

    print("Computation time: %.0f seconds." % (time.time() - start))

    #print(fit.x)
    #print(fit.lines)
    #print(fit.groups)
    
    # Plot results
    fig = plt.figure(figsize = (10, 8))
    fig.canvas.set_window_title('Lines')
    fig.suptitle(name)
    gs = gridspec.GridSpec(1, 1) 
    ax_0 = fig.add_subplot(gs[0])
    gs.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
    ax_0.plot(spec.x, spec.y)
    ax_0.plot(conv.x, conv.y)
    ax_0.scatter(lines.x, lines.y)

    plt.show()

if __name__ == '__main__':
    main()
