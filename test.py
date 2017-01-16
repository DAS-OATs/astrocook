
#Import libraries
from astrocook import spec1d
from astrocook import spec1dreader
from astropy import units as u
from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt
import copy

#prova

#Enable automatic unit display
from astropy.visualization import quantity_support
quantity_support()

#Create some fake data
a = [1,2,3]

#Create a spec1d object providing x, y, dy and units informations
s = spec1d(a, a, a, xUnit=u.Angstrom, yUnit=1.e-17 * u.erg / u.second / u.cm**2 / u.Angstrom)

#Test unit conversion capabilities
print()
print("Wavelength conversions:")
print(s.x * s.xUnit)
s.convert(xUnit=u.Hz)
print(s.x * s.xUnit)
s.convert(xUnit=u.Angstrom)
print(s.x * s.xUnit)

print()
print("Flux density conversions:")
print(s.y * s.yUnit)
s.convert(yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Hz)
print(s.y * s.yUnit)
s.convert(yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom)
print(s.y * s.yUnit)

print(s.y * s.yUnit)
s.convert(yUnit=u.watt / u.cm**2 / u.Angstrom)
print(s.y * s.yUnit)
s.convert(yUnit=1.e-17*u.erg / u.second / u.cm**2 / u.Angstrom)
print(s.y * s.yUnit)



#Plot the spec1d data
plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.show()



#Create a spec1dreader object
r = spec1dreader()

#Read data from a FITS file
s = r.sdss_dr10('spec-0752-52251-0323.fits')

#Plot the spec1d data
plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.title("MJD="+str(s.meta['MJD']) + ", plate="+str(s.meta['PLATEID']) + ", fiber="+str(s.meta['FIBERID']))
plt.show()


#Plot only the "good" spectral channels
s.useGood = True

plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.title("MJD="+str(s.meta['MJD']) + ", plate="+str(s.meta['PLATEID']) + ", fiber="+str(s.meta['FIBERID']))
plt.show()


#Create a powerlaw model
guess_pl = models.PowerLaw1D(amplitude=10, x_0=s.x.mean(), alpha=0.5)

#Plot data and model 
plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.plot(s.x * s.xUnit, guess_pl(s.x) * s.yUnit)
plt.show()

#Fit data
fit_t = fitting.LevMarLSQFitter()
fit_pl = fit_t(guess_pl, s.x, s.y)
print("Fit results:")
print(fit_pl)

#Plot data and model 
plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.plot(s.x * s.xUnit, fit_pl(s.x) * s.yUnit)
plt.show()


#De-redden spectrum
dered = copy.deepcopy(s)
dered.deredden(3.1 * 0.06846)

plt.plot(s.x * s.xUnit, s.y * s.yUnit)
plt.plot(dered.x * dered.xUnit, dered.y * dered.yUnit)
plt.show()
