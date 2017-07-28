from astrocook import *
from astropy import units as u
import matplotlib.pyplot as plt

xunit = u.nm
yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)

def main():

    # Input data
    name = 'J0940_CIV2'

    # Read the 1D spectrum
    spec = Spec1DReader().uves(name + '_spec.fits')

    # Create a Line object from the spectrum
    line = Line(spec)

    # Find lines in the spectrum
    line.find(kappa=5.0)

    # Plot lines
    line.plot()


if __name__ == '__main__':
    main()
