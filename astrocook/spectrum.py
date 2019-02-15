from .frame import Frame
from .line_list import LineList
from .message import *
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

    def convolve_gauss(self, std=20, input_col='y', output_col='conv'):
        """@brief Convolve a spectrum column with a profile using FFT transform.
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

    def extract_nodes(self, delta_x=1500, xunit=au.km/au.s):
        """ @brief Extract nodes from a spectrum. Nodes are averages of x and y
        in slices, computed after masking lines.
        @param delta_x Size of slices
        @param xunit Unit of wavelength or velocity
        @return Spectrum with nodes
        """
        try:
            xunit = au.Unit(xunit)
            delta_x = float(delta_x)*xunit
        except:
            print(prefix, msg_param_fail)

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

        self._nodes = Spectrum(x, xmin, xmax, y, dy, self._xunit, self._yunit)
        return self._nodes

    def interp_nodes(self, smooth=0):
        """ @brief Interpolate nodes with a univariate spline to estimate the
        emission level.
        @param smooth Smoothing of the spline
        @return 0
        """
        if not hasattr(self, '_nodes'):
            print(prefix, "I need nodes to interpolate. Please try Spectrum > "
                  "Extract nodes first.")
            return None

        try:
            smooth = float(smooth)
        except:
            print(prefix, msg_param_fail)

        x = self._nodes.x.value
        y = self._nodes.y.value
        dy = self._nodes.dy.value
        spl = uspline(x, y, w=dy, s=smooth)
        cont = spl(self.x)*self._yunit
        print(prefix, "I'm using interpolation as continuum.")
        if 'cont' in self._t.colnames:
            print(prefix, "I'm updating column 'cont'.")
        else:
            print(prefix, "I'm adding column 'cont'.")
        self._t['cont'] = cont #spl(self.x)
        self._lines._t['cont'] = np.interp(self._lines.x, self.x, cont)
        self._nodes._t['cont'] = np.interp(self._nodes.x, self.x, cont)
        return 0

    def find_peaks(self, col='conv', kind='min', kappa=3.0):
        """ @brief Find the peaks in a spectrum column. Peaks are the extrema
        (minima or maxima) that are more prominent than a given number of
        standard deviations. They are saved as a list of lines.
        @param col Column where to look for peaks
        @param kind Kind of extrema ('min' or 'max')
        @param kappa Number of standard deviations
        @return 0
        """

        if col not in self.t.colnames:
            print(prefix, "The spectrum has not a column named '%s'. Please "\
                  "pick another one." % col)
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
        lines = ext._copy(sel)
        self._lines = LineList(lines.x, lines.xmin, lines.xmax, lines.y,
                               lines.dy, self._xunit, self._yunit, self._meta)
        #self._lines = ext._copy(sel)
        self._lines_kind = 'peaks'
        self._mask_lines()
        print(prefix, "I'm using peaks as lines.")

        return self._lines
