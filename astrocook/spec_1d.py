from . import List
from .utils import many_gauss, savitzky_golay, dict_wave
from astropy import units as u
from astropy.constants import c
from astropy.io import fits as fits
from astropy.table import Column, Table
import copy
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve
from specutils import extinction

class Spec1D():
    """Class for generic spectra
    
    A generic spectrum is a Table with the following columns: 
        -# @x: channels;
        -# @xmin: lower limit for each channel; 
        -# @xmax: upper limit for each channel; 
        -# @y: flux density in the channel;
        -# @dy: error on @y.
        -# @group: quality/grouping of the channel;
        -# @resol: spectral resolution in the channel; 
        
    The following metadata are also expected:
        -# @order: spectrum order;
        -# @exptime: integration time of all exposures;
        -# @meta: miscellaneous information (TBD).        
    """
    
    def __init__(self, x, y, 
                 xmin=None,
                 xmax=None,
                 dx=None,
                 dy=None,
                 group=None,
                 resol=None,
                 xunit=None, #u.dimensionless_unscaled, 
                 yunit=None, #u.dimensionless_unscaled, 
                 exptime=float('nan'), 
                 order=-1,
                 meta=None,
                 dtype=float):
        ''' Constructor for the Spec1D class. '''

        col_x  = Column(np.asarray(copy.deepcopy(x) , dtype=float), name='X')
        col_y  = Column(np.asarray(copy.deepcopy(y) , dtype=float), name='Y')

        if (xmin is None   and   xmax is not None):
            raise Exception('When XMAX is used also XMIN must be given')
        if (xmin is not None   and   xmax is None):
            raise Exception('When XMIN is used also XMAX must be given')
        if (xmin is None   and   xmax is None):
            if (dx is None): 
                dx = np.roll((x - np.roll(x, 1))/2., -1)
                dx[-1] = dx[-2]
            xmin = x - dx
            xmax = x + dx
        col_xmin = Column(np.asarray(copy.deepcopy(xmin), dtype=dtype), name='XMIN')
        col_xmax = Column(np.asarray(copy.deepcopy(xmax), dtype=dtype), name='XMAX')

        if (dy is None):
            dy = np.repeat(float('nan'), len(col_x))
        col_dy = Column(np.asarray(copy.deepcopy(dy), dtype=dtype), name='DY')
        
        if (group is None):
            group = np.ones(len(col_x))
        col_g = Column(np.asarray(copy.deepcopy(group), dtype=int), name='GROUP')
        
        if (resol is None):
            resol = np.repeat(float('nan'), len(col_x))
        col_r = Column(np.asarray(copy.deepcopy(resol), dtype=dtype), name='RESOL')
        
        # Auxiliary data and meta
        self._exptime = float(exptime)
        self._order   = int(order)
        if (meta is None):
            meta = {}

        # Table creation
        self._t = Table(
            data=(col_xmin, col_xmax, col_x, col_y, col_dy, col_g, col_r), 
            masked=True, meta=meta)
        self._t['XMIN'].unit = xunit
        self._t['XMAX'].unit = xunit
        self._t['X'].unit    = xunit
        self._t['Y'].unit    = yunit
        self._t['DY'].unit   = yunit

        self._t['XMIN'].mask  = np.isnan(self._t['XMIN'].quantity.value)
        self._t['XMAX'].mask  = np.isnan(self._t['XMAX'].quantity.value)
        self._t['X'].mask     = np.isnan(self._t['X'].quantity.value)
        self._t['Y'].mask     = np.isnan(self._t['Y'].quantity.value)
        self._t['DY'].mask    = np.isnan(self._t['DY'].quantity.value)
        self._t['RESOL'].mask = np.isnan(self._t['RESOL'].quantity.value)

        self._useGood = False


    def _getWithMask(self,colName):
        if self._useGood:
            ret = self._t[colName].quantity[self._igood]
        else:
            ret = self._t[colName].quantity
        null = np.argwhere(self._t[colName].mask)
        if null.size > 0:
            ret[null] = float('nan')
        return ret

    @property
    def t(self):
        if self._useGood:
            return self._t[self._igood]
        else:
            return self._t

    @property
    def useGood(self):
        """Tells whether x, y, etc. getters return only data from channels flagged as good."""
        return self._useGood

    @useGood.setter
    def useGood(self, value):
        self._useGood = value
        if self._useGood:
            self._igood = np.argwhere(self._t['GROUP'] >= 0)

    @property
    def exptime(self):
        """Spectrum exposure time (in seconds)."""
        return self._exptime

    @property
    def order(self):
        """Spectrum order (-1 if not specified)."""
        return self._order

    @property
    def meta(self):
        """Meta information for the spectrum."""
        return self._t.meta

    def nchan(self):
        """Number of channels in the spectrum."""
        if self._useGood:
            return len(self._igood)
        else:
            return len(self._t)

    @property
    def x(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('X')

    @x.setter
    def x(self, value):
        if self._useGood:
            self._t['X'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['X'] = np.asarray(value, dtype='float')
        #self._t['X'].unit = self._t['XMIN'].unit

    @property
    def xmin(self):
        """Lower limit of quantity associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('XMIN')

    @xmin.setter
    def xmin(self, value):
        if self._useGood:
            self._t['XMIN'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['XMIN'] = np.asarray(value, dtype='float')
        #self._t['XMIN'].unit = self._t['X'].unit

    @property
    def xmax(self):
        """Upper limit of quantity associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self._getWithMask('XMAX')

    @xmax.setter
    def xmax(self, value):
        if self._useGood:
            self._t['XMAX'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['XMAX'] = np.asarray(value, dtype='float')
        #self._t['XMAX'].unit = self._t['X'].unit

    @property
    def dx(self):
        """Widths of spectrum channels, in the same units as the spectrum channels."""
        xmin = self._getWithMask('XMIN')
        xmax = self._getWithMask('XMAX')
        return xmax - xmin

    @property
    def y(self):
        """Quantities associated to spectrum intensities (e.g. flux density, luminosity density, nuF_nu, lambdaF_lambda, etc.)."""
        return self._getWithMask('Y')

    @y.setter
    def y(self, value):
        if self._useGood:
            self._t['Y'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['Y'] = np.asarray(value, dtype='float')
        #self._t['Y'].unit = self._t['DY'].unit

    @property
    def dy(self):
        """Uncertainties associated to y values."""
        return self._getWithMask('DY')

    @dy.setter
    def dy(self, value):
        if self._useGood:
            self._t['DY'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['DY'] = np.asarray(value, dtype='float')
        #self._t['DY'].unit = self._t['Y'].unit

    @property
    def group(self):
        """Return group flag for each spectrum channel."""
        if self._useGood:
            return self._t['GROUP'].quantity.value[self._igood]
        else:
            return self._t['GROUP'].quantity.value

    @group.setter
    def group(self, value):
        if self._useGood:
            self._t['GROUP'][self._iGood] = np.asarray(value, dtype='int')
        else:
            self._t['GROUP'] = np.asarray(value, dtype='int')

    @property
    def resol(self):
        """Return resolution for each spectrum channel."""
        return self._getWithMask('RESOL')

    @property
    def xunit(self):
        """Physical unit for the x property, to be expressed as an astropy unit."""
        return self._t['X'].unit

    @property
    def yunit(self):
        """Physical unit for the y property, to be expressed as an astropy unit."""
        return self._t['Y'].unit

    def convert(self, xunit=None, yunit=None):
        """Convert x and/or y values into equivalent quantities."""
        if not (xunit is None):
            mask = self._t['X'].mask
            q = self._t['X']
            p = q.to(xunit, equivalencies=u.spectral())
            self._t['X'] = p
            self._t['X'].mask = mask

            mask = self._t['XMIN'].mask
            q = self._t['XMIN']
            p = q.to(xunit, equivalencies=u.spectral())
            self._t['XMIN'] = p
            self._t['XMIN'].mask = mask

            mask = self._t['XMAX'].mask
            q = self._t['XMAX']
            p = q.to(xunit, equivalencies=u.spectral())
            self._t['XMAX'] = p
            self._t['XMAX'].mask = mask

        if not (yunit is None):
            mask = self._t['Y'].mask
            q = self._t['Y']
            p = q.to(yunit, equivalencies=u.spectral_density(self._t['X']))
            self._t['Y'] = p
            self._t['Y'].mask = mask

            mask = self._t['DY']
            q = self._t['DY']
            p = q.to(yunit, equivalencies=u.spectral_density(self._t['X']))
            self._t['DY'] = p
            self._t['DY'].mask = mask


    def todo_convert_logx(self, xunit=None, yunit=None):
        """Convert x and/or y values into equivalent quantities."""
        if not (xunit is None):
            #TODO: check following code
            if xunit.is_equivalent(u.m / u.s) \
                or self.xunit.is_equivalent(u.m / u.s):
                    equiv = [(u.nm, u.km / u.s, 
                              lambda x: np.log(x) * c.to(u.km / u.s).value, 
                              lambda x: np.exp(x / c.to(u.km / u.s).value)   \
                          )]

            mask = self._t['X'].mask
            q = self._t['X']
            p = q.to(xunit, equivalencies=equiv)
            self._t['X'] = p
            self._t['X'].mask = mask

            mask = self._t['XMIN'].mask
            q = self._t['XMIN']
            p = q.to(xunit, equivalencies=equiv)
            self._t['XMIN'] = p
            self._t['XMIN'].mask = mask

            mask = self._t['XMAX'].mask
            q = self._t['XMAX']
            p = q.to(xunit, equivalencies=equiv)
            self._t['XMAX'] = p
            self._t['XMAX'].mask = mask



    def convolve(self, col='y', prof=None, gauss_sigma=20, mode='same'):
        """Convolve a spectrum with a profile using FFT transform
        
        The profile must have the same length of the column X of the spectrum.
        If no profile @a prof is provided, a normalized Gaussian profile with 
        @a gauss_sigma is applied."""

        conv = copy.deepcopy(self)
        conv_col = getattr(conv, col)
        conv.todo_convert_logx(xunit=u.km/u.s)        
        if prof is None:
            par = np.stack([[np.mean(conv.x.value)], [gauss_sigma], [1]])
            par = np.ndarray.flatten(par, order='F')
            prof = many_gauss(conv.x.value, *par)
            prof = prof / np.sum(prof)
        conv.y = fftconvolve(conv_col, prof, mode=mode)
        conv.todo_convert_logx(xunit=self.x.unit)
        return conv
        
    def deredden(self, A_v, model='od94'):
        extFactor = extinction.reddening(self._t['X'], A_v, model=model)
        self._t['Y']  *= extFactor
        self._t['DY'] *= extFactor
    
    def find_extrema(self):
        """Find the extrema in a spectrum and save them as a spectrum"""
    
        min_idx = np.hstack(argrelmin(self.y))
        max_idx = np.hstack(argrelmax(self.y))
        extr_idx = np.sort(np.append(min_idx, max_idx))
        minima = self.from_table(self._t[min_idx])
        maxima = self.from_table(self._t[max_idx])
        extrema = self.from_table(self._t[extr_idx])
        return minima, maxima, extrema
        
    def find_lines(self, mode='abs', diff='max', kappa=3.0, hwidth=2):
        """Find the lines in a spectrum and save them as a line list"""
        
        # Find the extrema
        minima, maxima, extr = self.find_extrema()
        
        # Compute the difference between each extremum and its neighbours
        # N.B. Negative fluxes give wrong results! To be fixed 
        diff_y_left = (extr.y[:-2] - extr.y[1:-1]) 
        diff_y_right = (extr.y[2:] - extr.y[1:-1]) 
        if mode is 'em':
            diff_y_left = -diff_y_left
            diff_y_right = -diff_y_right
        
        # Check if the difference is above threshold
        diff_y_min = np.minimum(diff_y_left, diff_y_right)
        diff_y_max = np.maximum(diff_y_left, diff_y_right)
        thres = extr.dy[1:-1] * kappa
        if diff is 'min':
            line_pos = np.greater(diff_y_min, thres)
        else:
            line_pos = np.greater(diff_y_max, thres)
        line_x = extr.x[1:-1][line_pos]
        line_y = extr.y[1:-1][line_pos]

        if len(line_x) > 0: 

            # Compute the boundaries for line fitting    
            window = extr.rolling_window(hwidth * 2 + 1)
            bound = np.full(hwidth + 1, np.amin(self.x))
            if mode is 'em':
                bound = np.append(bound, 
                    window.x[np.arange(len(window.x)), np.argmin(window.y, 1)])
            else:
                bound = np.append(bound,
                    window.x[np.arange(len(window.x)), np.argmax(window.y, 1)])
            bound = np.append(bound, np.full(hwidth + 1, np.amax(self.x)))
            line_xmin = bound[:-hwidth * 2][line_pos]
            line_xmax = bound[hwidth * 2:][line_pos]
            
            # Create the linelist
            list = List(line_x, line_y, xmin=line_xmin, xmax=line_xmax)
            return list

    def from_table(self, table, meta = {}):
        """Read a spectrum from a (spectrum-like) table"""
        
        xmin = table['XMIN']
        xmax = table['XMAX']
        x = table['X']
        y = table['Y']            
        dy = table['DY']
        group = table['GROUP']
        resol = table['RESOL']
        dx = 0.5 * (xmax - xmin)
        
        c1 = np.argwhere(y > 0)
        c2 = np.argwhere(dy > 0)
        igood = np.intersect1d(c1, c2)

        good = np.repeat(-1, len(x))
        good[igood] = 1

        spec = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xunit=x.unit, 
                      yunit=y.unit, group=good, resol=resol, meta=meta)
        return spec
        
    def rolling_window(self, width):
        """Convert a spectrum in rolling-window format
    
        Rolling-window format means that each column entry is an array of size
        @p width, obtained by rolling a window through each column. The size
        of the spectrum is 
        
        Code adapted from 
        http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html"""

        shape = self.x.shape[:-1] + (self.x.shape[-1] - width + 1, width)
        strides = self.x.strides + (self.x.strides[-1],)
        window_xmin = np.lib.stride_tricks.as_strided(
            self.xmin, shape=shape, strides=strides)
        window_xmax = np.lib.stride_tricks.as_strided(
            self.xmax, shape=shape, strides=strides)
        window_x = np.lib.stride_tricks.as_strided(
            self.x, shape=shape, strides=strides)
        window_y = np.lib.stride_tricks.as_strided(
            self.y, shape=shape, strides=strides)
        window_dy = np.lib.stride_tricks.as_strided(
            self.dy, shape=shape, strides=strides)
        window_group = np.lib.stride_tricks.as_strided(
            self.group, shape=shape, strides=strides)
        window_resol = np.lib.stride_tricks.as_strided(
            self.resol, shape=shape, strides=strides)
        window = Spec1D(
            window_x, window_y, xmin=window_xmin, xmax=window_xmax,
            dy=window_dy, group=window_group, resol=window_resol)
        #window = Spec1DReader().table(window_t)
        return window
        
    def save(self, filename):
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', array=self.xmin),
             fits.Column(name='XMAX', format='E', array=self.xmax),
             fits.Column(name='X', format='E', array=self.x),
             fits.Column(name='Y', format='E', array=self.y),
             fits.Column(name='DY', format='E', array=self.dy),
             fits.Column(name='GROUP', format='I', array=self.group),
             fits.Column(name='RESOL', format='E', array=self.resol)])
        hdu.writeto(filename, overwrite=True)

    def to_wave(self, ion):
        """ Convert line redshifts to wavelengths using a reference """

        if (len(ion) == 1):
            ion = np.full(len(self.t), ion[0])
        
        if ((len(ion) > 1) and (len(ion) != len(self.t))):
            raise Exception("Ion list and line list must have the same length.")
            
        ion = np.asarray(ion)
        ion_ravel = np.ravel(ion, 'F')
        ion_wave = np.asarray([dict_wave[i].value for i in ion_ravel]) \
                   * dict_wave['Ly_a'].unit
        z = np.ravel(self.x.value, 'F')
        zmin = np.ravel(self.xmin.value, 'F') 
        zmax = np.ravel(self.xmax.value, 'F')
        x = (z + 1) * ion_wave
        xmin = (zmin + 1) * ion_wave
        xmax = (zmax + 1) * ion_wave
        self.x = np.unique(x)
        self.xmin = np.unique(xmin)
        self.xmax = np.unique(xmax)

    def to_z(self, ion):
        """ Convert line wavelengths to redshifts using a reference """

        if (len(ion) == 1):
            ion = np.full(len(self.t), ion[0])

        if ((len(ion) > 1) and (len(ion) != len(self.t))):
            raise Exception("Ion list and line list must have the same length.")
            
        ion = np.asarray(ion)
        ion_ravel = np.ravel(ion, 'F')
        ion_wave = np.asarray([dict_wave[i].value for i in ion_ravel]) \
                   * dict_wave['Ly_a'].unit
        x = np.resize(self.x.value, ion_ravel.shape) * self.x.unit
        xmin = np.resize(self.xmin.value, ion_ravel.shape) * self.xmin.unit
        xmax = np.resize(self.xmax.value, ion_ravel.shape) * self.xmax.unit
        z = x / ion_wave - 1
        zmin = xmin / ion_wave - 1
        zmax = xmax / ion_wave - 1
        self.x = np.reshape(z, ion.shape, 'F')
        self.xmin = np.reshape(zmin, ion.shape, 'F')
        self.xmax = np.reshape(zmax, ion.shape, 'F')

