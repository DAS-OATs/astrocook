from .frame import Frame
from .linelist import LineList
from .message import *
from astropy import units as au
#from astropy import constants as aconst
#from astropy import table as at
from copy import deepcopy as dc
#from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve
from scipy.interpolate import UnivariateSpline as uspline

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

    def _mask_lines(self):
        """ @brief Create a mask consisting on the ['xmin', 'xmax'] regions from
        the associated line list
        @return 0
        """

        x = self._safe(self.x)
        mask = np.zeros(len(x), dtype=bool)
        for (xmin, xmax) in zip(self._lines.xmin, self._lines.xmax):
            #xmin = l['xmin']*self._lines._xunit
            #xmax = l['xmin']*self._lines._xunit
            mask += np.logical_and(x>xmin, x<xmax)
        if 'lines_mask' in self._t.colnames:
            print(prefix, "I'm updating column 'lines_mask'.")
        else:
            print(prefix, "I'm adding column 'lines_mask'.")
            self._t['lines_mask'] = np.empty(len(self.x), dtype=bool)
        self._t['lines_mask'][self._where_safe] = mask
        
        return 0

    def convolve_gauss(self, std=20, input_col='y', output_col='conv'):
        """@brief Convolve a spectrum colum  with a profile using FFT transform.
        @param std Standard deviation of the gaussian (km/s)
        @param input_col Input column
        @param output_col Output column
        @return 0
        """

        try:
            std = float(std) * au.km/au.s
        except:
            print(prefix, msg_param_fail)

        # Create profile
        xunit = self.x.unit
        self._convert_x()
        x = self._safe(self.x)
        mean = np.median(x)
        prof = np.exp(-((x - mean) / std).value**2)
        if (len(prof) % 2 == 0):
            prof = prof[:-1]
        prof = prof / np.sum(prof)

        # Convolve
        if output_col in self._t.colnames:
            print(prefix, "I'm updating column '%s'." % output_col)
        else:
            print(prefix, "I'm adding column '%s'." % output_col)
        conv = dc(self._t[input_col])
        safe = self._safe(conv)
        conv[self._where_safe] = fftconvolve(safe, prof, mode='same')\
                                              *self._t[input_col].unit
        self._t[output_col] = conv
        self._convert_x(xunit=xunit)

        return 0

    def interp_emission(self, s=0):
        """ @brief Estimate the emission level by interpolating across lines.
        @param s Spline smoothing
        @return 0
        """

        s = float(s)

        if not hasattr(self, '_lines'):
            print(prefix, "I need lines to interpolate the emission. Please "
                  "try Spectrum > Find peaks.")
            return None
        self._xmask(self._lines)
        """
        x = self.x[self._t['xmask']].value
        y = self.y[self._t['xmask']].value
        dy = self.dy[self._t['xmask']].value
        spl = uspline(x, y, w=dy, s=s)
        self._t['cont'] = spl(self.x)
        if 'cont' in self._t.colnames:
            print(prefix, "I'm updating column 'cont'.")
        else:
            print(prefix, "I'm adding column 'cont'.")
        print(prefix, "I'm using interpolation as continuum.")
        """
        return 0

    def find_peaks(self, col='conv', kind='min', kappa=3.0):
        """ @brief Find the peak in a spectrum column. Peaks are the extrema
        (minima or maxima) that are more prominent than a given number of
        standard deviations. They are saved as a list of lines.
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
            ext.xmin = np.append(self.x[0], ext.x[:-1])
            ext.xmax = np.append(ext.x[1:], self.x[-1])

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
        self._lines = LineList()
        self._lines = ext._copy(sel)
        self._lines_kind = 'peaks'
        self._mask_lines()
        print(prefix, "I'm using peaks as lines.")

        return self._lines
