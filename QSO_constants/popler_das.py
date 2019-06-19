from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

gc = fits.open('J1103-2645_RSPEC_PRE_spec.fits')
mm = fits.open('j1103m2645_0.7_error2_spec.fits')
gct = Table(gc[1].data)
mmt = Table(mm[1].data)

plt.plot(gct['x'], gct['dy']/gct['cont'], label='DAS error')
plt.plot(mmt['x']/10, mmt['dy'], label='Popler error')
plt.yscale('log')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux  (normalized)')
plt.legend()
plt.show()


dy_int = np.interp(gct['x'], mmt['x']/10, mmt['dy'])
plt.plot(gct['x'], dy_int/(gct['dy']/gct['cont']), label='Popler/DAS error')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux  (ratio)')
plt.legend()
plt.show()

def window_rms_2(a, window_size):
    ret = np.zeros(len(a))
    for i in range(window_size//2-1, len(a)-window_size//2):
        asel = a[i-window_size//2+1: i+window_size//2]
        ret[i] = np.std(asel)
    return ret

gc_rms = window_rms_2(gct['y'], 50)
mm_rms = window_rms_2(mmt['y'], 50)

plt.plot(gct['x'], gc_rms/gct['cont'], label='DAS RMS', color='C2')
plt.plot(mmt['x']/10, mm_rms, label='Popler RMS', color='C3')
plt.plot(gct['x'], gct['dy']/gct['cont'], label='DAS error', color='C0')
plt.plot(mmt['x']/10, mmt['dy'], label='Popler error', color='C1')
plt.yscale('log')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux  (normalized)')
plt.legend()
plt.show()

plt.plot(gct['x'], gc_rms/gct['cont'], label='DAS RMS', color='C2')
plt.plot(gct['x'], gct['dy']/gct['cont'], label='DAS error', color='C0')
plt.yscale('log')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux  (normalized)')
plt.legend()
plt.show()

plt.plot(mmt['x']/10, mm_rms, label='Popler RMS', color='C3')
plt.plot(mmt['x']/10, mmt['dy'], label='Popler error (expected fluctuation)', color='C1')
plt.yscale('log')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux  (normalized)')
plt.legend()
plt.show()
