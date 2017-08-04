from astrocook import *
from astropy import units as u
import numpy as np
import sys

yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def main():

    # Input data
    #name = 'J0940_CIV_1'  # Single CIV system 
    #name = 'J0940_CIV_2'  # Another single CIV system
    name = 'J0940_CIV_forest'  # Whole chunk of CIV forest (time consuming)

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Find lines in the spectrum
    line = Line(spec)
    line.find(kappa=5.0)
    
    # Create a Syst object from the lines
    syst = Syst(line, spec, ion='CIV')

    # Create redshift table
    syst.create_z()
    
    # Match redshifts
    syst.match_z()

    # Flatten redshifts
    syst.flatten_z()
    
    ltot = len(syst.t)
    for l in range(ltot):
        print("Line %i of %i (%3.2f)..." % (l + 1, ltot, syst.x[l].value),
              end=" ", flush=True)
        
        # Define the line group around a given redshift
        group = syst.group(line=l)
        
        # Define the spectral chunk around a given redshift
        chunk = syst.chunk(line=l)

        # Guess the unabsorbed continuum
        unabs_guess = syst.unabs(group, chunk)
    
        # Guess the Voigt profile
        voigt_guess = syst.voigt(group, chunk)

        # Fit the model 
        fit = syst.fit(group, chunk, unabs_guess, voigt_guess)

        # Fit the model, incrementally adding components
        # NOT WORKING YET!
        #group, chunk = syst.auto(line=l)
    
        # Plot lines and system
        print("close graph to continue.")
        #syst.plot(group, chunk)  # To visually match components
        syst.plot(group, chunk, split=True)  # To inspect the fit  
    
if __name__ == '__main__':
    main()
