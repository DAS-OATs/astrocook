from astrocook import *
from astropy import units as u
from copy import deepcopy as dc
import matplotlib.pyplot as plt
import numpy as np

xunit = u.nm
yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def main():

    filename = 'B2126-15_Lya_f_8_spec.fits'
    filename = 'B1122_all_spec.fits'
    spec = Spec1DReader().uves(filename)
    #spec.plot()
    line = Line(spec)
    line.find(kappa=5.0, sigma=20.0)
    #line.plot2()
    line.cont(smooth=1.0)
    #line.plot2()
    ltot = len(line.t)
    line_i = dc(line)
    x_arr = line_i.x
    line._resol = 60000
    group_check = 0
    for l in range(ltot):
        print("Line %i of %i (%3.2f %s)..." \
              % (l + 1, ltot, x_arr[l].value, x_arr[l].unit),
              end=" ", flush=True)
        if (np.array_equal(line_i.group(x=x_arr[l])[1], group_check)):
            print("same group, skipping.")
            pass
        else:
            group_check = line_i.group(x=x_arr[l])[1]
            #group = line.group(x=x_arr[l])
            #chunk = line.chunk(x=x_arr[l])
            #line.plot2(group, chunk)
            #psf = line.psf(group, chunk, line._resol)
            #norm_guess = line.norm(group, chunk)
            #voigt_guess = line.voigt(group, chunk)
            #line.plot2(group, chunk)
            #fit = line.fit(group, chunk, norm_guess, voigt_guess, psf)
            #line.plot2(group, chunk)
            group, chunk = line.auto(x=x_arr[l])
            #line.plot2(group, chunk)
    line.plot2()            
        
if __name__ == '__main__':
    main()
