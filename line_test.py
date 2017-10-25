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
    #name = 'J0940_Lya_f_2'    # Whole chunk of Lya forest (time consuming)
    #name = 'B2126-15_Lya_f_7' # Another whole chunk of Lya forest
    #name = 'B1937_all' # A whole spectrum    
    
    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')
    
    # Create a Line object from the spectrum
    line = Line(spec)
    line._resol = 60000
    
    # Find absorption lines in the spectrum
    line.find(kappa=5.0, sigma=20.0)

    # Estimate the continuum
    line.cont()


    line.fit_list(plot=True)
    
    """
    ltot = len(line.t)
    line_i = dc(line)
    x_arr = line_i.x
    group_check = 0
    #for l in range(14,15):
    for l in range(ltot):
        
        print("Line %i of %i (%3.2f %s)..." \
              % (l + 1, ltot, x_arr[l].value, x_arr[l].unit),
              end=" ", flush=True)
    
        # Check if the group is new
        if (np.array_equal(line_i.group(x=x_arr[l])[1], group_check)):
            print("same group, skipping.")
            pass
        else:
            group_check = line_i.group(x=x_arr[l])[1]

            # Fit the model, incrementally adding components
            print(l)
            group, chunk = line.fit_auto(x=x_arr[l], i_max=1)
    
            # Plot lines
            print("Close graph to continue.")
            line.plot(group, chunk)

    #line.plot()
    """

if __name__ == '__main__':
    main()
