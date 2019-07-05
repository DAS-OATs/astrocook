from .frame import Frame
from .line_list import LineList
#from .syst_list import SystList
from .message import *
#from .vars import *
from astropy import units as au
#from astropy import constants as aconst
#from astropy import table as at
from copy import deepcopy as dc
#from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve
from scipy.interpolate import UnivariateSpline as uspline
from scipy.stats import sem

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
                 dtype=float,
                 cont=[]):
        super(Spectrum, self).__init__(x, xmin, xmax, y, dy, xunit, yunit, meta,
                                       dtype)
        if cont != []:
            self._t['cont'] = cont*self._yunit

    def _copy(self, sel=None):
        copy = super(Spectrum, self)._copy(sel)
        cols = [c for c in self._t.colnames \
                if c not in ['x', 'xmin', 'xmax', 'y', 'dy']]
        for c in cols:
            copy._t[c] = self._t[c][sel]
        return copy

    def _convolve_gauss(self, std=20, input_col='y', output_col='conv',
                        verb=True):

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
        if verb:
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

    def _extract_nodes(self, delta_x=1500, xunit=au.km/au.s):

        self._slice(delta_x, xunit)

        x_ave = []
        xmin_ave = []
        xmax_ave = []
        y_ave = []
        dy_ave = []
        if 'lines_mask' not in self._t.colnames:
            print(prefix, "Lines weren't masked. I'm taking all spectrum.")

        for s in self._slice_range:
            try:
                where_s = np.where(np.logical_and(self._t['slice']==s,
                                                  self._t['lines_mask']==0))
            except:
                where_s = np.where(self._t['slice']==s)

            #print(self.x[np.where(self._t['slice']==s)][0], len(where_s[0]))
            if len(where_s[0])>0:
                x_where_s = self.x[where_s].value
                y_where_s = self.y[where_s].value
                dy_where_s = self.dy[where_s].value
                x_ave.append(np.average(x_where_s))
                xmin_ave.append(x_where_s[0])
                xmax_ave.append(x_where_s[-1])
                y_ave.append(np.average(y_where_s, weights=dy_where_s))
                dy_ave.append(sem(y_where_s))
        x = np.array(x_ave) * self._xunit
        xmin = np.array(xmin_ave) * self._xunit
        xmax = np.array(xmax_ave) * self._xunit
        y = np.array(y_ave) * self._yunit
        dy = np.array(dy_ave) * self._yunit

        return x, xmin, xmax, y, dy

    def _interp_nodes(self, lines, nodes, smooth=0):
        """ @brief Interpolate nodes with a univariate spline to estimate the
        emission level.
        @param smooth Smoothing of the spline
        @return 0
        """

        x = nodes.x.value
        y = nodes.y.value
        dy = nodes.dy.value
        spl = uspline(x, y, w=dy, s=smooth)
        cont = spl(self.x)*self._yunit
        print(prefix, "I'm using interpolation as continuum.")
        if 'cont' in self._t.colnames:
            print(prefix, "I'm updating column 'cont'.")
        else:
            print(prefix, "I'm adding column 'cont'.")
        self._t['cont'] = cont #spl(self.x)
        lines._t['cont'] = np.interp(lines.x, self.x, cont)
        nodes._t['cont'] = np.interp(nodes.x, self.x, cont)
        return 0

    def _find_peaks(self, col='conv', kind='min', kappa=3.0, **kwargs):

        y = self._safe(self._t[col])
        min_idx = np.hstack(argrelmin(y, **kwargs))
        max_idx = np.hstack(argrelmax(y, **kwargs))
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
        lines = ext._copy(sel)

        return lines

    def fit_syst(self, series='CIV', z=1.6971, logN=13, b=10, resol=45000):
        """ @brief Add a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @return System list
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        self._systs = SystList(self)
        self._systs._add_fit(series, z, logN, b, resol)

        return self._systs

    def _mask_lines(self, lines):
        """ @brief Create a mask consisting on the ['xmin', 'xmax'] regions from
        the associated line list
        @return 0
        """

        x = self._safe(self.x)
        mask = np.zeros(len(x), dtype=bool)
        for (xmin, xmax) in zip(lines.xmin, lines.xmax):
            mask += np.logical_and(x>=xmin, x<=xmax)
        if 'lines_mask' in self._t.colnames:
            print(prefix, "I'm updating column 'lines_mask'.")
        else:
            print(prefix, "I'm adding column 'lines_mask'.")
            self._t['lines_mask'] = np.empty(len(self.x), dtype=bool)
        self._t['lines_mask'][self._where_safe] = mask

        return 0

    def _slice(self, delta_x=1000, xunit=au.km/au.s):
        """ @brief Create 'slice' columns. 'slice' columns contains an
        increasing counter to split 'x' values into evenly-sized slices
        (typically defined in velocity space).
        @param delta_x Size of slices
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        xunit_orig = self._xunit
        self._convert_x(xunit=xunit)
        x = self._safe(self.x)
        self._t['slice'] = np.empty(len(self.x), dtype=int)
        self._t['slice'][self._where_safe] = np.array(x//delta_x)
        self._slice_range = range(self._t['slice'][self._where_safe][0],
                                  self._t['slice'][self._where_safe][-1])
        self._convert_x(xunit=xunit_orig)
        return 0
