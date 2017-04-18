#Import libraries
from astrocook import Spec1D
from astrocook import Spec1DReader
from astropy import units as u
from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt
import copy



#Enable automatic unit display
from astropy.visualization import quantity_support
quantity_support()

#Create some fake data
a = [1,2,3]

#Create a Spec1D object providing x, y, dy and units informations
s = Spec1D(a, a, dy=a, xUnit=u.Angstrom, yUnit=1.e-17 * u.erg / u.second / u.cm**2 / u.Angstrom)

#Test unit conversion capabilities
print()
print("Wavelength conversions:")
print(s.x)
s.convert(xUnit=u.Hz)
print(s.x)
s.convert(xUnit=u.Angstrom)
print(s.x)

print()
print("Flux density conversions:")
print(s.y)
s.convert(yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Hz)
print(s.y)
s.convert(yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom)
print(s.y)

print(s.y)
s.convert(yUnit=u.watt / u.cm**2 / u.Angstrom)
print(s.y)
s.convert(yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom)
print(s.y)


#Plot the Spec1D data
plt.plot(s.x, s.y)
plt.show()



#Create a Spec1DReader object
r = Spec1DReader()

#Read data from a FITS file
s = r.sdss_dr10('spec-0752-52251-0323.fits')

#Plot the Spec1D data
plt.plot(s.x, s.y)
plt.title("MJD="+str(s.meta['MJD']) + ", plate="+str(s.meta['PLATEID']) + ", fiber="+str(s.meta['FIBERID']))
plt.show()


#Plot only the "good" spectral channels
s.useGood = True

plt.plot(s.x, s.y)
plt.title("MJD="+str(s.meta['MJD']) + ", plate="+str(s.meta['PLATEID']) + ", fiber="+str(s.meta['FIBERID']))
plt.show()


#Create a powerlaw model
guess_pl = models.PowerLaw1D(amplitude=10, x_0=s.x.mean(), alpha=0.5)

#Plot data and model 
plt.plot(s.x, s.y)
plt.plot(s.x, guess_pl(s.x).value * s.yUnit)
plt.show()

#Fit data
fit_t = fitting.LevMarLSQFitter()
fit_pl = fit_t(guess_pl, s.x.value, s.y.value)
print("Fit results:")
print(fit_pl)

#Plot data and model 
plt.plot(s.x, s.y)
plt.plot(s.x, fit_pl(s.x.value) * s.yUnit)
plt.show()


#De-redden spectrum
dered = copy.deepcopy(s)
dered.deredden(3.1 * 0.06846)

plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.plot(dered.x * dered.xUnit, dered.y * dered.yUnit)
plt.show()
