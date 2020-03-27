from .frame import Frame
from .line_list import LineList
#from .syst_list import SystList
from .message import *
#from .vars import *
from astropy import units as au
from astropy.modeling.models import BlackBody
from astropy.modeling.powerlaws import PowerLaw1D
#from astropy import constants as aconst
#from astropy import table as at
from copy import deepcopy as dc
import logging
#from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve
from scipy.interpolate import UnivariateSpline as uspline
from scipy.stats import sem
from tqdm import tqdm

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
                 cont=[],
                 resol=[]):
        super(Spectrum, self).__init__(x, xmin, xmax, y, dy, xunit, yunit, meta,
                                       dtype)
        if cont != []:
            self._t['cont'] = cont*self._yunit
        if resol != []:
            self._t['resol'] = resol


    def _copy(self, sel=None):
        copy = super(Spectrum, self)._copy(sel)
        cols = [c for c in self._t.colnames \
                if c not in ['x', 'xmin', 'xmax', 'y', 'dy']]
        for c in cols:
            copy._t[c] = self._t[c][sel]
        return copy

    def _gauss_convolve(self, std=20, input_col='y', output_col='conv',
                        verb=True):

        # Create profile
        xunit = self.x.unit
        self._x_convert()
        x = self._safe(self.x)
        mean = np.median(x)
        prof = np.exp(-((x - mean) / std).value**2)
        if (len(prof) % 2 == 0):
            prof = prof[:-1]
        prof = prof / np.sum(prof)

        # Convolve
        if verb:
            if output_col not in self._t.colnames:
                logging.info("I'm adding column '%s'." % output_col)
        conv = dc(self._t[input_col])
        safe = self._safe(conv)
        conv[self._where_safe] = fftconvolve(safe, prof, mode='same')\
                                              *self._t[input_col].unit
        self._t[output_col] = conv
        self._x_convert(xunit=xunit)

        return 0


    """
    def _syst_fit(self, series='CIV', z=1.6971, logN=13, b=10, resol=45000):
        @brief Add a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @return System list

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        self._systs = SystList(self)
        self._systs._add_fit(series, z, logN, b, resol)

        return self._systs
    """


    def _lines_mask(self, lines, source=None):
        """ @brief Create a mask consisting on the ['xmin', 'xmax'] regions from
        the associated line list
        @return 0
        """

        x = self._safe(self.x)
        mask = np.zeros(len(x), dtype=bool)

        if source is not None:
            where = lines.t['source'] == source
        else:
            where = range(len(lines.t))

        lines_xmin = lines.xmin[where]
        lines_xmax = lines.xmax[where]

        for (xmin, xmax) in zip(lines_xmin, lines_xmax):
            mask += np.logical_and(x>=xmin, x<=xmax)
        if 'lines_mask' not in self._t.colnames:
            logging.info("I'm adding column 'lines_mask' to spectrum.")
            self._t['lines_mask'] = np.empty(len(self.x), dtype=bool)
        self._t['lines_mask'][self._where_safe] = mask

        return 0


    def _nodes_clean(self, nodes, kappa=5.0):
        """
        y_sel = [True]
        hw = window//2
        while np.any(y_sel):
            y_start = nodes.y.value
            yg = [y_start[max(0, i-hw):i+hw+1] for i in range(len(y_start))]
            yg_median = [np.median(g) for g in yg]
            yg_std = [np.std(g) for g in yg]

            y_sel = [np.abs(y-m) > kappa*s \
                     for y, m, s in zip(y_start, yg_median, yg_std)]
            #print(yg)
            #print(yg_median)
            #print(yg_std)
            #print(y_sel)
            print(np.sum(y_sel))
            #print(nodes._t[y_sel])
            nodes._t.remove_rows(y_sel)
        """
        lines = nodes._peaks_find(col='y', kappa=kappa, kind='max')
        y_sel = np.where(np.in1d(nodes.x.value, lines.x.value))[0]
        nodes._t.remove_rows(y_sel)
        while np.any(y_sel):
            lines = nodes._peaks_find(col='y', kappa=kappa)
            y_sel = np.where(np.in1d(nodes.x.value, lines.x.value))[0]
            #print(len(y_sel))
            nodes._t.remove_rows(y_sel)
        #"""
        return nodes


    def _nodes_extract(self, delta_x=1500, xunit=au.km/au.s):

        self._slice(delta_x, xunit)
        x_ave = []
        xmin_ave = []
        xmax_ave = []
        y_ave = []
        dy_ave = []
        if 'lines_mask' not in self._t.colnames:
            logging.warning("Lines weren't masked. I'm taking all spectrum.")

        for s in self._slice_range:
            try:
                where_s = np.where(np.logical_and(self._t['slice']==s,
                                                  self._t['lines_mask']==0))
            except:
                where_s = np.where(self._t['slice']==s)

            # Use deabs column if present
            y = self._t['deabs'] if 'deabs' in self._t.colnames else self.y.value

            if len(where_s[0])>0.1*len(np.where(self._t['slice']==s)[0]):
                x_where_s = self.x[where_s].value
                y_where_s = y[where_s]
                dy_where_s = self.dy[where_s].value
                x_ave.append(np.median(x_where_s))
                xmin_ave.append(x_where_s[0])
                xmax_ave.append(x_where_s[-1])
                y_ave.append(np.median(y_where_s))
                dy_ave.append(sem(y_where_s))
        x = np.array(x_ave) * self._xunit
        xmin = np.array(xmin_ave) * self._xunit
        xmax = np.array(xmax_ave) * self._xunit
        y = np.array(y_ave) * self._yunit
        dy = np.array(dy_ave) * self._yunit

        return Spectrum(x, xmin, xmax, y, dy, self._xunit, self._yunit)


    def _nodes_interp(self, lines, nodes, smooth=0):
        """ @brief Interpolate nodes with a univariate spline to estimate the
        emission level.
        @param smooth Smoothing of the spline
        @return 0
        """


        x = nodes.x.value
        y = nodes.y.value
        dy = nodes.dy.value
        isnan = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),
                              np.isnan(dy))
        spl = uspline(x[~isnan], y[~isnan], w=dy[~isnan], s=smooth)
        cont = spl(self.x)*self._yunit
        logging.info("I'm using interpolation as continuum.")
        if 'cont' not in self._t.colnames:
            logging.info("I'm adding column 'cont'.")
        self._t['cont'] = cont #spl(self.x)
        lines._t['cont'] = np.interp(lines.x, self.x, cont)
        nodes._t['cont'] = np.interp(nodes.x, self.x, cont)
        return 0


    def _peaks_find(self, col='conv', kind='min', kappa=3.0, **kwargs):

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
        #for m,M,l,r in zip(ext.xmin, ext.xmax, diff_y_left, diff_y_right):

        #    print(m,M,l,r)
        diff_y_max = np.minimum(diff_y_left, diff_y_right)

        # +1 is needed because sel is referred to the [1:-1] range of rows
        # in the spectrum
        sel = np.where(np.greater(diff_y_max, ext.dy[1:-1] * kappa))[0]+1
        lines = ext._copy(sel)

        return lines

    def _rebin(self, dx, xunit, y, dy):

        # Convert spectrum into chosen unit
        # A deep copy is created, so the original spectrum is preserved

        self.t.sort('x')
        self._x_convert(xunit=xunit)

        # Create x, xmin, and xmax
        from .format import Format
        format = Format()
        xstart, xend = np.nanmin(self.x), np.nanmax(self.x)
        x = np.arange(xstart.value, xend.value, dx) * xunit
        xmin, xmax = format._create_xmin_xmax(x)

        # Compute y and dy combining contributions
        im = 0
        iM = 1
        xmin_in = self.xmin[iM].value
        xmax_in = self.xmax[im].value
        y_out = np.array([])
        dy_out = np.array([])
        for i, (m, M) \
            in enum_tqdm(zip(xmin.value, xmax.value), len(xmin),
                         "spectrum: Rebinning"):
            while xmin_in < M:
                iM += 1
                try:
                    xmin_in = self.xmin[iM].value
                except:
                    break
            while xmax_in < m:
                im += 1
                try:
                    xmax_in = self.xmax[im].value
                except:
                    break
            ysel = y[im:iM+1]
            dysel = dy[im:iM+1]
            y_out = np.append(y_out, np.average(ysel, weights=1/dysel**2))
            dy_out = np.append(dy_out, np.sqrt(np.sum(dysel**2/dysel**4))\
                                               /np.sum(1/dysel**2))

        # Create a new spectrum and convert it to the units of the original one
        out = Spectrum(x, xmin, xmax, y_out, dy_out, xunit=xunit, yunit=y.unit,
                       meta=self.meta)
        out._x_convert(xunit=self._xunit_old)
        return out


    def _resol_est(self, px, update):
        dx = self.xmax[px-1:].value-self.xmin[:1-px].value
        dx = np.append(np.append([dx[0]]*(px//2), dx), [dx[-1]*(px//2)])\
             *self.x.unit
        if update:
            if 'resol' not in self._t.colnames:
                logging.info("I'm adding column 'resol'.")
            else:
                logging.info("I'm updating column 'resol'.")
            self._t['resol'] = self.x/dx
        logging.info('The mean estimated resolution is %4.2f.'
                     % np.mean(self._t['resol']))
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
        self._x_convert(xunit=xunit)
        x = self._safe(self.x)
        self._t['slice'] = np.empty(len(self.x), dtype=int)
        self._t['slice'][self._where_safe] = np.array(x//delta_x)
        self._slice_range = range(self._t['slice'][self._where_safe][0],
                                  self._t['slice'][self._where_safe][-1])
        self._x_convert(xunit=xunit_orig)
        return 0

    def _template_bb(self, temp=6000, scale=1.0):
        bb = BlackBody(temperature=temp*au.K, scale=scale*au.erg/(self._xunit*au.cm**2*au.s*au.sr))
        output_col = 'blackbody'
        if output_col not in self._t.colnames:
            logging.info("I'm adding column '%s'." % output_col)
        self._t[output_col] = bb(self.x)


    def _template_pl(self, ampl=1.0, x_ref=None, index=-1.0):
        if x_ref == None: x_ref = np.mean(self.x)
        pl = PowerLaw1D(amplitude=ampl, x_0=x_ref, alpha=-index)
        output_col = 'power_law'
        if output_col not in self._t.colnames:
            logging.info("I'm adding column '%s'." % output_col)
        self._t[output_col] = pl(self.x)



    def _zap(self, xmin, xmax):

        w = np.where(np.logical_and(self.x.value>xmin, self.x.value<xmax))
        self._t['y'][w] = np.interp(
                              self.x[w].value,
                              [self.x[w][0].value, self.x[w][-1].value],
                              [self.y[w][0].value, self.y[w][-1].value])*self._yunit
        return 0
