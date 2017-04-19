from astrocook import List, ListReader, Spec1D, Spec1DCont, Spec1DReader
from astrocook.utils import savitzky_golay
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import copy
from cycler import cycler
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import sys
import time

import warnings
warnings.filterwarnings("ignore")

def find_lines(spec, col='y', sigma=10, mode='abs', kappa=4.0, diff='max'):
    """Specialized version of Spec1D.find_lines, with prior convolution
    
    The convolution algorithm is provisional, that's why it's not been moved to
    Spec1D."""
    
    """
    print("  Rebinning spectrum...")
    rebin_x = np.trunc(spec.x / bin)
    grouped = spec.t.group_by(rebin_x)
    binned = grouped.groups.aggregate(np.mean)
    rebin = spec.from_table(binned)
    rebin_copy = copy.deepcopy(rebin)
    attr = getattr(rebin_copy, col)
    rebin_copy.y = attr
    """
    print("  Filtering spectrum...")
    conv = spec.convolve(col=col, gauss_sigma=sigma)

    print("  Finding lines...")
    #list = rebin_copy.find_lines(mode=mode, kappa=kappa, diff=diff)
    list = conv.find_lines(mode=mode, kappa=kappa, diff=diff)
    if list is not None:
        print("  %i lines found." % (len(list.t)))
    else:
        print("  No lines found.")
    return list, conv

def main():
    """Determine the continuum of a QSO spectrum"""

    start = time.time()
    timer = 10.0
    
    name = str(sys.argv[1])
    
    # Set up plotting window
    fig = plt.figure(figsize = (12, 8))
    fig.canvas.set_window_title('Continuum')
    fig.suptitle(name)
    gs = gridspec.GridSpec(5, 1) 
    ax_0 = fig.add_subplot(gs[0])
    ax_1 = fig.add_subplot(gs[1])
    ax_2 = fig.add_subplot(gs[2])
    ax_3 = fig.add_subplot(gs[3])   
    ax_4 = fig.add_subplot(gs[4])   
    ax_0.set_prop_cycle(cycler('color', ['blue', 'red', 'cyan']))
    ax_1.set_prop_cycle(cycler('color', ['blue', 'red', 'cyan']))
    ax_2.set_prop_cycle(cycler('color', ['blue', 'red', 'green']))
    ax_3.set_prop_cycle(cycler('color', ['blue', 'green']))
    ax_4.set_prop_cycle(cycler('color', ['blue', 'red', 'green', 'cyan', 'black']))
    ax_3.set_ylim([0.9, 1.1])   
    gs.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
    
    spec_read = Spec1DReader()
    list_read = ListReader()

    try:
        spec = spec_read.uves(name + '_spec.fits')
    except:
        spec = spec_read.simul(name + '_spec.fits')
    print("Original spectrum loaded.")
        
        
    # Compute the smoothing windows
    # N.B. The input spectrum should have a fixed velocity bin
    smooth = 500 * u.km / u.s
    spec.todo_convert_logx(xUnit=u.km/u.s)
    narrow = int(smooth / np.mean(spec.xmax - spec.xmin))
    medium = narrow * 2
    wide = narrow * 4
    narrow = narrow + 1 if narrow % 2 == 0 else narrow
    medium = medium + 1 if medium % 2 == 0 else medium
    wide = wide + 1 if wide % 2 == 0 else wide
    spec.todo_convert_logx(xUnit=u.nm)    
    print("Savitzky-Golay smoothing windows: ")
    print(" wide (absorption-clean continuum): %i pixels" % (wide))
    print(" narrow (emission continuum): %i pixels" % (narrow))
    print(" medium (final continuum): %i pixels" % (medium))
    
    # Remove emission lines
    try: 

        spec_abs = spec_read.cont(name + '_spec_abs.fits')
        list_abs = list_read.gen(name + '_list_abs.fits')
        spec_em = spec_read.cont(name + '_spec_em.fits')
        list_em = list_read.gen(name + '_list_em.fits')
        print("Clean spectrum loaded.")

    except:
    
        # Remove absorption lines
        try: 

            spec_abs = spec_read.cont(name + '_spec_abs.fits')
            list_abs = list_read.gen(name + '_list_abs.fits')
            print("Absorption-clean spectrum loaded.")
        
        except:
            
            print("Clean spectrum from absorption lines:")
            try:
                list_abs = list_read.gen(name + '_list_abs.fits')
                print(" Absorption lines loaded.")
            except:
                print(" Extracting absorption lines...")
                list_abs, conv_abs = find_lines(spec, kappa=10.0)
                if (list_abs is not None):
                    list_abs.save(name + '_list_abs.fits')

            print(" Grouping absorption lines...")
            list_abs.group()

            print(" Fitting and removing absorption lines...")
            spec_abs = copy.deepcopy(spec)
            #spec_abs.save(name + '_spec_abs.fits')
            spec_abs = Spec1DCont(spec)
            #spec_abs.conv = copy.deepcopy(conv_abs.y)
            #spec_abs.fit_lines(list_abs, method='sg', sg_size=51)
            spec_abs.fit_lines(list_abs, timer=0.1)

            print("Fitting a Savitzky-Golay continuum to the spectrum...")
            spec_abs.sg(spec_abs.abs_rem, wide, 3)
            spec_abs.save(name + '_spec_abs.fits')
            
        #ax_2.plot(spec_abs.x, spec_abs.y / spec_abs.cont, 'r-')
    
        print("Clean spectrum from emission lines:")
        try:
            list_em = list_read.gen(name + '_list_em.fits')
            print(" Emission lines loaded.")
        except:
            print(" Extracting emission lines...")
            list_em, conv_em = find_lines(spec_abs, col='abs_rem', sigma=300,
                                           mode='em', kappa=10.0, diff='min')
            if list_em is not None:
                list_em.save(name + '_list_em.fits')

        if list_em is not None:
            print(" Grouping emission lines...")
            list_em.group()

        print(" Fitting and removing emission lines...")
        spec_em = copy.deepcopy(spec_abs)
        spec_em.fit_lines(list_em, mode='em', col='abs_rem', method='sg', 
                          sg_size=narrow, timer=0.1)
        #spec_em.conv = copy.deepcopy(conv_abs.y)
        spec_em.cont = copy.deepcopy(spec_em.em_fit)
        spec_em.save(name + '_spec_em.fits')


    try:
        spec_pre = spec_read.cont(name + '_spec_pre.fits')
        print("Preliminary continuum loaded.")
    except:
        print("Saving the preliminary continuum...")
        # The continuum is extracted by taking the em_fit profile where the emission
        # lines are and filtering the abs_rem profile elsewhere
        spec_pre = Spec1DCont(spec)
        spec_pre.abs_fit = copy.deepcopy(spec_abs.abs_fit)
        spec_pre.abs_rem = copy.deepcopy(spec_abs.abs_rem)
        spec_pre.em_fit = copy.deepcopy(spec_em.em_fit)
        spec_pre.em_rem = copy.deepcopy(spec_em.em_rem)
        #spec_pre_conv = spec_pre.convolve(col='abs_rem', gauss_sigma=1000)
        #spec_pre.cont = copy.deepcopy(spec_pre_conv.y)
        spec_pre.cont = copy.deepcopy(spec_abs.cont)
        where_em = np.where(~np.isnan(spec_pre.em_fit))
        spec_pre.cont[where_em] = spec_pre.em_fit[where_em]
        spec_pre.save(name + '_spec_pre.fits')
    
    try:
        spec_cont = spec_read.cont(name + '_spec_cont.fits')
        spec_res = spec_read.cont(name + '_spec_res.fits')
        list_cont = list_read.gen(name + '_list_cont.fits')
        print("Final continuum loaded.")
    except:
        print("Extract the final continuum...")
        spec_cont = copy.deepcopy(spec_pre)
        try:
            list_cont = list_read.gen(name + '_list_cont.fits')
            print(" Lines to be finally removed loaded.")
        except:
            print(" Extracting lines to be finally removed...")
            list_cont, conv_cont = find_lines(spec_cont, kappa=5.0, sigma=5.0)
            list_cont.save(name + '_list_cont.fits')

        print(" Grouping absorption lines...")
        list_cont.group()

        print(" Fitting and removing absorption lines...")
        spec_cont.fit_lines(list_cont, method='v', timer=timer)
        
        print("Fitting a Savitzky-Golay continuum to the spectrum...")
        spec_cont.sg(spec_cont.abs_rem, medium, 3)

        """
        print(" Extracting residuals...")
        spec_res = copy.deepcopy(spec_cont)
        spec_res.y = spec_res.abs_rem
        list_res, conv_res = find_lines(spec_res, kappa=2.0, sigma=5.0)
        print(" Grouping residuals...")
        list_res.group()

        print(" Fitting residuals...")
        spec_res.fit_lines(list_res, method='v', timer=timer)

        print("Fitting a Savitzky-Golay continuum to the spectrum...")
        spec_res.sg(spec_res.abs_rem, medium, 3)
        spec_cont.save(name + '_spec_res.fits')

        spec_cont.cont = spec_res.cont
        """

        spec_cont.save(name + '_spec_cont.fits')

    print("Computation time: %.0f seconds." % (time.time() - start))
    
    # Plot results
    ax_0.plot(spec.x, spec.y) #blue
    ax_0.plot(spec_abs.x, spec_abs.abs_fit) #red
    #ax_0.plot(spec_abs.x, spec_abs.conv) #cyan
    ax_0.scatter(list_abs.x, list_abs.y, c='black', marker='+')
    ax_0.set_ylabel("Y")
    ax_1.plot(spec_abs.x, spec_abs.abs_rem) #blue
    ax_1.plot(spec_abs.x, spec_abs.cont) #red
    ax_1.plot(spec_em.x, spec_em.em_fit) #cyan
    ax_1.set_ylabel("Y")
    #ax_1.plot(rebin_em.x, rebin_em.y)
    if list_em is not None:
        ax_1.scatter(list_em.x, list_em.y, s=400, c='black', marker='+')
    ax_2.plot(spec.x, spec.y) #blue
    ax_2.plot(spec_pre.x, spec_pre.cont) #red
    ax_2.plot(spec_cont.x, spec_cont.cont) #green
    ax_2.set_ylabel("Y")
    ax_3.plot(spec.x, spec.y / spec_cont.cont) #blue
    ax_3.plot(spec_cont.x, spec_cont.cont / spec_cont.cont) #green
    ax_3.set_ylabel("Y / CONT")
#    ax_3.set_xlabel("X (" + "{}".format(spec.x.unit) + ")")
    ax_4.plot(spec_cont.x, spec_cont.y) #blue
    ax_4.scatter(list_cont.x, list_cont.y, c='black', marker='+')
    ax_4.plot(spec_abs.x, spec_abs.abs_rem) #red
    ax_4.plot(spec_cont.x, spec_cont.abs_fit) #green 
    ax_4.plot(spec_cont.x, spec_cont.abs_rem) #cyan
    #ax_4.plot(spec_res.x, spec_res.abs_rem) #black
    ax_4.set_ylabel("Y / CONT")
    ax_4.set_xlabel("X (" + "{}".format(spec.x.unit) + ")")
    plt.show()

    fig = plt.figure(figsize = (12, 10))
    plt.plot(spec_cont.x, spec_cont.y) #blue
    plt.scatter(list_cont.x, list_cont.y, c='black', marker='+')
    plt.plot(spec_abs.x, spec_abs.abs_rem) #red
    plt.plot(spec_cont.x, spec_cont.abs_fit) #green 
    plt.plot(spec_cont.x, spec_cont.abs_rem) #cyan
    #ax_4.plot(spec_res.x, spec_res.abs_rem) #black
    plt.show()

    """
    # Ad usum Stephani, for the H2020 slide 
    fig_show = plt.figure(figsize = (15, 5))
    gs = gridspec.GridSpec(1, 1) 
    ax_0 = fig_show.add_subplot(gs[0])
    ax_0.set_prop_cycle(cycler('color', ['black', 'green', 'blue']))
    #gs.tight_layout(fig_show, rect=[0.01, 0.01, 1, 0.97])
    ax_0.plot(spec.x, spec.y)
    ax_0.plot(spec_cont.x, spec_cont.abs_fit)
    ax_0.scatter(list_abs.x, list_abs.y, c='red', marker='+')
    ax_0.plot(spec_cont.x, spec_cont.cont)
    ax_0.set_ylabel("Flux")
    ax_0.set_xlabel("Wavelength (" + "{}".format(spec.x.unit) + ")")
    plt.show()
    """
    
    print("Success!")

if __name__ == '__main__':
    main()
