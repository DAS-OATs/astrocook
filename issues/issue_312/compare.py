from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

base = Table(fits.open('A_best_systs.fits')[1].data)
rfit = Table(fits.open('A_best_refit_systs.fits')[1].data)
newc = Table(fits.open('A_best_refit_newconv_systs.fits')[1].data)

base.sort(['z0', 'series'])
rfit.sort(['z0', 'series'])
newc.sort(['z0', 'series'])

x = rfit
y = newc

col = 'z'
dcol = 'dz'

#plt.errorbar(x[col], y[col], xerr=x[dcol], yerr=y[dcol], fmt='o', alpha=0.2)
#plt.plot([-1000,1000], [-1000,1000], color='black')
#plt.ylim(np.min(y[col]), np.max(y[col]))

plt.scatter(x[col], y[col]-x[col], alpha=0.2)
plt.plot([-1000,1000], [0, 0], color='black')

plt.xlim(np.min(x[col]), np.max(x[col]))
plt.show()