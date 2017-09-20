from astrocook import *
from astropy import units as u
from copy import deepcopy as dc
import numpy as np

def main():

    plot = True
    name = 'J0940_Lya_f_2'
    #name = 'B2126-15_Lya_f_2' 
    
    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Create a Line object from the spectrum
    line = Line(spec)
    line._resol = 60000

    # Find absorption lines in the spectrum
    line.find(kappa=5.0, sigma=20.0)

    # Estimate the continuum
    line.cont()

    # Create a model of the lines and guess the column density
    model = Model(spec, line=line)
    N_guess = model.N_guess(line._cont)
    N_argsort = np.argsort(N_guess)
    #print(N_guess[N_argsort])

    ltot = len(line.t)
    #for l in range(ltot):
    N_fit = np.array([])

    flux_corr = line.corr_tau(N_thres=1e30 / u.cm**2)
    spec_rem = dc(spec)
    for i in range(3):#ltot):

        N_argsel = N_argsort[-i-2:-1]
        pl = 's'                
        if (i == 0):
            pl = ''
        print("%i line group%s out of %i..." \
              % (i + 1, pl, ltot), end=" ", flush=True)
    
        spec_rem.cont(flux_corr=flux_corr)
        line_i = dc(line)
        line_i._cont = dc(spec_rem._cont)

        x_arr = line_i.x
        for x in x_arr[N_argsel]:
            group = line_i.group(x)
            chunk = line_i.chunk(x)
            psf = line_i.psf(group, chunk, line_i._resol)
            norm_guess = line_i.norm(group, chunk)
            voigt_guess = line_i.voigt(group, chunk)
            cont_guess = norm_guess
            fit = line_i.fit(group, chunk, cont_guess, voigt_guess, psf)
            N_fit = np.append(N_fit, line_i._N_fit)

        spec_rem.y = line_i._rem.y
        N_thres = np.median(N_fit)
        flux_corr = line_i.corr_tau(N_thres)
        print(N_thres, flux_corr)
                   
        if (plot == True):
            print("Close graphs to continue.")
            line_i.plot(block=False)
            spec_rem.plot(block=False)
            line_i.plot_N(N_fit)
        else:
            print("")

    line_i.plot(block=False)
    line_i.plot_N(N_fit)            
if __name__ == '__main__':
    main()
