from astrocook import *
from astropy import units as u
from copy import deepcopy as dc
import numpy as np
import sys

yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def main():

    # Input data
    name = 'J0940_CIV_1'  # Single CIV system 
    #name = 'J0940_CIV_2'  # Another single CIV system
    #name = 'J0940_CIV_f'  # Whole chunk of CIV forest (time consuming)

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Find lines in the spectrum
    line = Line(spec)
    line.find(kappa=5.0)
    
    # Create a Syst object from the lines
    syst = Syst(line, spec, ion='CIV')
    syst._resol = 60000
    
    # Estimate the continuum
    syst.cont(wind=10.0, low=1.5, fact=20)

    # Create redshift table
    syst.create_z()
    
    # Match redshifts
    syst.match_z()

    # Flatten redshifts
    syst.flatten_z()
    
    ltot = len(syst.t)
    syst_i = dc(syst)
    x_arr = syst_i.x
    group_check = 0
    for l in range(1):#ltot):
        print("Redshift %i of %i (%3.4f)..." % (l + 1, ltot, x_arr[l].value),
              end=" ", flush=True)

        # Check if the group is new
        if (np.array_equal(syst_i.group(x=x_arr[l])[1], group_check)):
            print("same group, skipping.")
        else:
            group_check = syst_i.group(x=x_arr[l])[1]
            
            #""" This part is commented because it is run automatically later
            # Define the line group around a given redshift
            group = syst.group(x=x_arr[l])
            
            # Define the spectral chunk around a given redshift
            chunk = syst.chunk(x=x_arr[l])

            # Model the instrumental PSF
            psf = line.psf(group, chunk, syst._resol)

            # First way: estimate continuum from scratch
            #unabs_guess = syst.unabs(group, chunk)
            #voigt_guess = syst.voigt(group, chunk)
            #cont_guess = unabs_guess

            # Second way: use existing continuum
            norm_guess = syst.norm(group, chunk)
            voigt_guess = syst.voigt(group, chunk)
            cont_guess = norm_guess
            
            # Fit the model 
            #fit = syst.fit(group, chunk, cont_guess, voigt_guess, psf)
            #"""

            """ This runs all the part above automatically """
            # Fit the model, incrementally adding components
            #group, chunk = syst.auto(x=x_arr[l])
    
            # Plot lines and system
            print("close graph to continue.")
            #syst.plot(group, chunk)  # To visually match components
            syst.plot(group, chunk, split=True)  # To inspect the fit

            

if __name__ == '__main__':
    main()
