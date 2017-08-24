from astrocook import *
from astropy import units as u
from copy import deepcopy as dc
import matplotlib.pyplot as plt
import numpy as np

xunit = u.nm
yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def main():

    # Input data
    name = 'J0940_Lya_1'    # Single Lya line complex
    #name = 'J0940_Lya_2'    # Another single Lya line complex
    #name = 'J0940_CIV_1'    # CIV system blindly fitted as Ly_a
    #name = 'J0940_Lya_f'    # Whole chunk of Lya forest (time consuming &
                             # not very effective because of continuum issues)
    #name = 'B2126-15_Lya_f' # Same as above
    
    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Create a Line object from the spectrum
    line = Line(spec)
    line._resol = 60000
    
    # Find absorption lines in the spectrum
    line.find(kappa=5.0)

    ltot = len(line.t)
    line_i = dc(line)
    x_arr = line_i.x
    group_check = 0
    for l in range(ltot):
        print("Line %i of %i (%3.2f %s)..." \
              % (l + 1, ltot, x_arr[l].value, x_arr[l].unit),
              end=" ", flush=True)
    
        # Check if the group is new
        if (np.array_equal(line_i.group(x=x_arr[l])[1], group_check)):
            print("same group, skipping.")
        else:
            group_check = line_i.group(x=x_arr[l])[1]

            """ This part is commented because it is run automatically later
            # Define the line group around a given wavelength
            group = line.group(x=x_arr[l])
            
            # Define the spectral chunk around a given wavelength
            chunk = line.chunk(x=x_arr[l])

            # Guess the unabsorbed continuum
            unabs_guess = line.unabs(group, chunk)

            # Guess the Voigt profile
            voigt_guess = line.voigt(group, chunk)

            # Model the instrumental PSF
            psf = line.psf(group, chunk, line._resol)

            # Fit the model 
            fit = line.fit(group, chunk, unabs_guess, voigt_guess, psf)
            """

            """ This runs all the part above automatically """
            # Fit the model, incrementally adding components
            group, chunk = line.auto(x=x_arr[l])
    
            # Plot lines
            print("close graph to continue.")
            line.plot(group, chunk)


if __name__ == '__main__':
    main()
