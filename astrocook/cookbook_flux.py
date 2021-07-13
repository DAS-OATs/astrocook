from .vars import filt_x_skymap, zero_point_skymap
from astropy import units as au
from astropy.io import fits
from copy import deepcopy as dc
import logging
from matplotlib import pyplot as plt
import numpy as np
import os

class CookbookFlux(object):
    """ Cookbook of utilities for flux calibration
    """

    def __init__(self):
        super(CookbookFlux, self).__init__()

        p = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) \
            + '/../filtskymap/'
        self._filt_name = p+'filt%sskymap.fits'
        self._filt_x = filt_x_skymap
        self._filt_xunit = au.Angstrom
        self._filt_yunit = 1e17*au.erg/au.s/au.cm**2/au.Angstrom
        self._zero_point = zero_point_skymap


    def _mags_compute(self, bands):
        """ @brief Compute magnitudes
        @details Compute the magnitudes of the spectrum in the SkyMapper bands.
        N.B. Zero-point flux in u band is ~6.87e-06 erg/s/cm2/A.
        @param bands List of bands ('u', 'v', 'g', 'r', 'i', or 'z')
        @return Magnitudes
        """

        mags = {}
        for band in bands:

            spec = dc(self.sess.spec)

            f = fits.open(self._filt_name % band)
            hdr = f[0].header
            data = f[0].data
            crval1 = hdr['CRVAL1']
            cdelt1 = hdr['CDELT1']
            naxis1 = hdr['NAXIS1']
            y = data * self._filt_yunit
            x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]

            xmin, xmax = np.min(x), np.max(x)
            spec_xmin, spec_xmax = np.min(spec.x.to(self._filt_xunit).value), \
                                   np.max(spec.x.to(self._filt_xunit).value)
            sel = np.where(np.logical_and(x>spec_xmin, x<spec_xmax))[0]
            xsel = x[sel]
            ysel = y[sel]
            if spec_xmin>xmin and spec_xmin<xmax:
                logging.warning("The red end of the %s filter falls outside the "
                "spectrum." % band)
            elif spec_xmax<xmax and spec_xmax>xmin:
                logging.warning("The blue end of the %s filter falls outside the "
                "spectrum." % band)
            elif len(sel)==0:
                logging.error("The %s filter does not overlap the spectrum."
                              % band)
            try:
                spec_r  = spec._rebin(xsel[0]*self._filt_xunit/spec.x.unit,
                                     (xsel[-1]+cdelt1)*self._filt_xunit/spec.x.unit,
                                     cdelt1, self._filt_xunit, spec.y, spec.dy)
                corr = np.sum(ysel)/np.sum(y)
                #plt.plot(spec_r.x.to(self._filt_xunit).value, spec_r.y*ysel/corr)
                #plt.plot(spec_r.x.to(self._filt_xunit).value, spec_r.y)
                #plt.plot(spec.x, spec.y)
                #plt.plot(xsel, ysel.value)
                #print("%3.8e" % np.median(spec_r.y.value))
                mag = -2.5 * np.log10(np.sum(spec_r.y.value*ysel.value)/corr) + self._zero_point[band]
                logging.info("AB magnitude in the %s filter: %3.2f."
                             % (band, mag))
                mags[band] = mag
            except:
                pass
        #plt.show()

        return mags


    def mags_adjust(self, bands, refs, deg=1):
        """ @brief Adjust magnitudes
        @details Adjust a spectrum to some reference magnitudes in the SkyMapper
        bands. The differences with the reference magnitudes are modeled with a
        polynomial regression to adjust the flux.
        @param bands List of bands ('u', 'v', 'g', 'r', 'i', or 'z')
        @param refs List of reference magnitudes
        @param deg Degree of polynomial regression
        @return Magnitudes
        """

        mags_in = self._mags_compute(bands)
        x = np.array([self._filt_x[b] for b in bands])
        diffs = np.array([float(refs[b])-mags_in[b] for b in bands])

        x = x[~np.isnan(diffs)]
        diffs = diffs[~np.isnan(diffs)]

        spec = self.sess.spec
        corr = 10**(np.poly1d(np.polyfit(x, diffs, deg)/-2.5)\
                (spec.x.to(au.nm).value))
        spec._t['y'] = spec._t['y']*corr
        spec._t['dy'] = spec._t['dy']*corr
        mags_in = self._mags_compute(bands)
