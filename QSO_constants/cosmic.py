from astropy.io import fits
import astroscrappy as ascr
import matplotlib.pyplot as plt
import numpy as np

test = fits.open('ESPRE.2019-05-23T00:17:14.397.fits')
data = test[1].data
med = np.median(data)
std = np.std(data)
mask, clean = ascr.detect_cosmics(data)

sel = np.s_[6200:7000,800:2200]


plt.imshow(data[sel], vmin=med-std, vmax=med+std)
#plt.imshow(clean[sel], vmin=med-std, vmax=med+std)
#plt.imshow(mask[sel])
plt.show()
