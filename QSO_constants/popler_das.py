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

mm = fits.open('j1103m2645_0.4.fits')
mm_h = mm[0].header
mm_wavelog = np.arange(mm_h['CRVAL1'], mm_h['CRVAL1']+mm_h['NAXIS1']*mm_h['CDELT1'], mm_h['CDELT1'])[:mm_h['NAXIS1']]
mm_wave = 10**mm_wavelog
mm_flux = mm[0].data[0]
mm_error = mm[0].data[1]
mm_expfluc = mm[0].data[2]
mm_rms = window_rms_2(mm_flux, 50)

plt.plot(mm_wave, mm_rms, label='RMS (51 pix)', color='C1')
plt.plot(mm_wave, mm_error, label='Error (row 2)', color='C0')
plt.plot(mm_wave, mm_expfluc, label='Exp. fluc. (row 3)', color='C3')
plt.legend()
plt.xlim((4750, 5000))
plt.ylim((0, 0.040))
plt.title(r'$\Delta$v = 0.4 km/s')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Normalised flux')
plt.show()
