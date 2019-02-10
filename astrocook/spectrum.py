from .frame import Frame
from .message import *
from astropy import units as au
from astropy import constants as aconst
from astropy import table as at
from copy import deepcopy as dc
#from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve

prefix = "Spectrum:"

class Spectrum(Frame):
    """Class for spectra

    A Spectrum is a Frame with methods for handling spectral operations."""

    def __init__(self,
                 x=[],
                 xmin=[],
                 xmax=[],
                 y=[],
                 dy=[],
                 xunit=au.nm,
                 yunit=au.erg/au.cm**2/au.s/au.nm,
                 meta={},
                 dtype=float):
        super(Spectrum, self).__init__(x, xmin, xmax, y, dy, xunit, yunit, meta,
                                       dtype)

    def _copy(self, sel=None):
        copy = super(Spectrum, self)._copy(sel)
        cols = [c for c in self._t.colnames \
                if c not in ['x', 'xmin', 'xmax', 'y', 'dy']]
        for c in cols:
            copy._t[c] = self._t[c][sel]
        return copy

    def convert_x(self, zem=0, xunit=au.km/au.s):
        """@brief Convert wavelengths into velocities and vice versa.
        @param zem Emission redshift
        @param xunit Unit of velocity or wavelength
        """

        try:
            zem = float(zem)
        except:
            print(prefix, msg_param_fail)
        xunit = au.Unit(xunit)

        xem = (1+zem) * 121.567*au.nm
        equiv = [(au.nm, au.km / au.s,
                  lambda x: np.log(x/xem.value)*aconst.c.to(au.km/au.s),
                  lambda x: np.exp(x/aconst.c.to(au.km/au.s).value)*xem.value)]

        self._xunit = xunit
        self.x = self.x.to(xunit, equivalencies=equiv)
        self.xmin = self.xmin.to(xunit, equivalencies=equiv)
        self.xmax = self.xmax.to(xunit, equivalencies=equiv)
        return 0

    def convolve_gauss(self, std=20):
        """@brief Convolve a spectrum with a profile using FFT transform. The
        convolution is saved in column 'conv'.
        @param std Standard deviation of the gaussian
        @return 0
        """

        try:
            std = float(std) * au.km/au.s
        except:
            print(prefix, msg_param_fail)

        # Create profile
        xunit = self.x.unit
        self.convert_x()
        x = self._safe(self.x)
        mean = np.median(x)
        prof = np.exp(-((x - mean) / std).value**2)
        if (len(prof) % 2 == 0):
            prof = prof[:-1]
        prof = prof / np.sum(prof)

        # Convolve
        if 'conv' not in self._t.colnames:
            print(prefix, "I'm adding column 'conv'.")
        yconv = dc(self.y)
        ysafe = self._safe(yconv)
        yconv[self._where_safe] = fftconvolve(ysafe, prof, mode='same')\
                                              *self.y.unit
        self._t['conv'] = yconv
        self.convert_x(xunit=xunit)

        return 0

    def extract_region(self, xmin, xmax):
        """ @brief Extract a spectral region as a new spectrum.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return Spectral region
        """

        try:
            xmin = float(xmin) * au.nm
            xmax = float(xmax) * au.nm
        except:
            print(prefix, msg_param_fail)
            return None

        if xmin > xmax:
            temp = xmin
            xmin = xmax
            xmax = temp
            print(prefix, msg_param_swap)


        reg = dc(self)
        reg.x.unit.to(au.nm)
        where = np.full(len(reg.x), True)
        s = np.where(np.logical_and(self._safe(reg.x) > xmin,
                                    self._safe(reg.x) < xmax))
        where[s] = False
        reg._t.remove_rows(where)

        if len(reg.t) == 0:
            print(prefix, msg_output_fail)
            return None
        else:
            return reg

    def find_peaks(self, col='conv', kind='min', kappa=3.0):
        """ @brief Find the peak in a spectrum column. Peaks are the extrema
        (minima or maxima) that are more prominent than a given number of
        standard deviations.
        @param col Column where to look for peaks
        @param kind Kind of extrema ('min' or 'max')
        @param kappa Number of standard deviations
        @return 0
        """

        if col not in self.t.colnames:
            print(prefix, "The spectrum has not a column named '%s'. Please "\
                  "pick another one" % col)
            return None
        kappa = float(kappa)

        y = self._safe(self._t[col])
        min_idx = np.hstack(argrelmin(y))
        max_idx = np.hstack(argrelmax(y))
        ext_idx = np.sort(np.append(min_idx, max_idx))
        ext = self._copy(ext_idx)

        if len(ext.t) > 0:
            # Set xmin and xmax from adjacent extrema
            ext.xmin[0] = self.x[0]
            ext.xmin[1:] = ext.x[:-1]
            ext.xmax[-1] = self.x[-1]
            ext.xmax[:-1] = ext.x[1:]

            diff_y_left = (ext._t[col][:-2] - ext._t[col][1:-1])
            diff_y_right = (ext._t[col][2:] - ext._t[col][1:-1])
            if kind == 'max':
                diff_y_left = -diff_y_left
                diff_y_right = -diff_y_right

        # Check if the difference is above threshold
        diff_y_max = np.maximum(diff_y_left, diff_y_right)
        # +1 is needed because sel is referred to the [1:-1] range of rows
        # in the spectrum
        sel = np.where(np.greater(diff_y_max, ext.dy[1:-1] * kappa))[0]+1
        self._peaks = ext._copy(sel)

        return 0
