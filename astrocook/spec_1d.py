from . import List
from .utils import convolve, convolve2, dict_doubl, dict_wave, dict_f, many_gauss, redchi_f, savitzky_golay
from astropy import units as u
from astropy.constants import c
from astropy.io import fits as fits
from astropy.table import Column, Table
import copy
from copy import deepcopy as dc
from lmfit import CompositeModel as lmc
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
    NavigationToolbar2TkAgg
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.signal import argrelmin, argrelmax, fftconvolve
from scipy.stats import sigmaclip
from specutils import extinction
#import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess
import time
#import tkinter as tk

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
                 xmin=[],
                 xmax=[],
                 dx=[],
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


        if ((xmin == []) != (xmax == [])):
            raise Exception("XMIN and XMAX must be provided together.")
        if ((xmin == []) and (xmax == [])):
            #if (np.size(x) > 0):
            if (x != []):
                if (dx == []):
                    dx = np.roll((x - np.roll(x, 1))/2., -1)
                    dx[-1] = dx[-2]
                xmin = x - dx
                xmax = x + dx

        data = ()
        if (x is not None):
            col_x  = Column(np.asarray(dc(x) , dtype=float), name='X')
            col_y  = Column(np.asarray(dc(y) , dtype=float), name='Y')
            data = (col_x, col_y)
            
        if (xmin != []):
            col_xmin = Column(np.asarray(dc(xmin), dtype=dtype), name='XMIN')
            col_xmax = Column(np.asarray(dc(xmax), dtype=dtype), name='XMAX')
            data += (col_xmin, col_xmax) 

        if data is ():
            data = None


        if (x is not None):
            if (dy is None):
                dy = np.repeat(float('nan'), len(col_x))
            col_dy = Column(np.asarray(dc(dy), dtype=dtype), name='DY')
            if (data is not None):
                data += (col_dy,)
            
            if (group is None):
                group = np.ones(len(col_x))
            col_group = Column(np.asarray(dc(group), dtype=int), name='GROUP')
            if (data is not None):
                data += (col_group,)
            
            if (resol is None):
                resol = np.repeat(float('nan'), len(col_x))
            col_resol = Column(np.asarray(dc(resol), dtype=dtype), name='RESOL')
            if (data is not None):
                data += (col_resol,)

        # Auxiliary data and meta
        self._exptime = float(exptime)
        self._order   = int(order)
        if (meta is None):
            meta = {}

        # Table creation
        self._t = Table(data=data, masked=True, meta=meta)
        if (x != []):
            self._t['X'].unit = xunit
            self._t['Y'].unit = yunit
            self._t['DY'].unit = yunit
            self._t['XMIN'].unit = xunit
            self._t['XMAX'].unit = xunit

            self._t['XMIN'].mask  = np.isnan(self._t['XMIN'].quantity.value)
            self._t['XMAX'].mask  = np.isnan(self._t['XMAX'].quantity.value)
            self._t['X'].mask     = np.isnan(self._t['X'].quantity.value)
            self._t['Y'].mask     = np.isnan(self._t['Y'].quantity.value)
            self._t['DY'].mask    = np.isnan(self._t['DY'].quantity.value)
            self._t['RESOL'].mask = np.isnan(self._t['RESOL'].quantity.value)

        self._use_good = False


    def _getWithMask(self,colName):
        if self._use_good:
            ret = self._t[colName].quantity[self._igood]
        else:
            ret = self._t[colName].quantity
        null = np.argwhere(self._t[colName].mask)
        if null.size > 0:
            ret[null] = float('nan')
        return ret

    @property
    def t(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t

    @property
    def use_good(self):
        return self._use_good

    @use_good.setter
    def use_good(self, value):
        self._use_good = value
        if self._use_good:
            self._igood = np.argwhere(self._t['GROUP'] >= 0)

    """Tells whether x, y, etc. getters return only data from channels flagged as good."""
    """
    @property
    def useGood(self):
        return self._useGood

    @useGood.setter
    def useGood(self, value):
        self._useGood = value
        if self._useGood:
            self._igood = np.argwhere(self._t['GROUP'] >= 0)
    """

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
        if self._use_good:
            return len(self._igood)
        else:
            return len(self._t)

    @property
    def x(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self.mask_col('X')

    @x.setter
    def x(self, value):
        if self._use_good:
            self._t['X'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['X'] = np.asarray(value, dtype='float')
        #self._t['X'].unit = self._t['XMIN'].unit

    @property
    def xmin(self):
        """Lower limit of quantity associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self.mask_col('XMIN')

    @xmin.setter
    def xmin(self, value):
        if self._use_good:
            self._t['XMIN'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['XMIN'] = np.asarray(value, dtype='float')
        #self._t['XMIN'].unit = self._t['X'].unit

    @property
    def xmax(self):
        """Upper limit of quantity associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        return self.mask_col('XMAX')

    @xmax.setter
    def xmax(self, value):
        if self._use_good:
            self._t['XMAX'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['XMAX'] = np.asarray(value, dtype='float')
        #self._t['XMAX'].unit = self._t['X'].unit

    @property
    def dx(self):
        """Widths of spectrum channels, in the same units as the spectrum channels."""
        xmin = self.mask_col('XMIN')
        xmax = self.mask_col('XMAX')
        return xmax - xmin

    @property
    def y(self):
        """Quantities associated to spectrum intensities (e.g. flux density, luminosity density, nuF_nu, lambdaF_lambda, etc.)."""
        return self.mask_col('Y')

    @y.setter
    def y(self, value):
        if self._use_good:
            self._t['Y'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['Y'] = np.asarray(value, dtype='float')
        #self._t['Y'].unit = self._t['DY'].unit

    @property
    def dy(self):
        """Uncertainties associated to y values."""
        return self.mask_col('DY')

    @dy.setter
    def dy(self, value):
        if self._use_good:
            self._t['DY'][self._iGood] = np.asarray(value, dtype='float')
        else:
            self._t['DY'] = np.asarray(value, dtype='float')
        #self._t['DY'].unit = self._t['Y'].unit

    @property
    def group(self):
        """Return group flag for each spectrum channel."""
        if self._use_good:
            return self._t['GROUP'].quantity.value[self._igood]
        else:
            return self._t['GROUP'].quantity.value

    @group.setter
    def group(self, value):
        if self._use_good:
            self._t['GROUP'][self._iGood] = np.asarray(value, dtype='int')
        else:
            self._t['GROUP'] = np.asarray(value, dtype='int')

    @property
    def resol(self):
        """Return resolution for each spectrum channel."""
        return self.mask_col('RESOL')

    @property
    def xunit(self):
        """Physical unit for the x property, to be expressed as an astropy unit."""
        return self._t['X'].unit

    @property
    def yunit(self):
        """Physical unit for the y property, to be expressed as an astropy unit."""
        return self._t['Y'].unit
    
    def cont(self, smooth=20.0, flux_corr=1.0):
        self._precont = dc(self)
        range_x = np.max(self.x) - np.min(self.x)
        x = self.x
        y = self.y * flux_corr
        self._cont = dc(self._precont)        
        clip_x = x
        clip_y = y
        stop = False
        i = 0
        frac = smooth*u.nm/range_x
        #le = sm.nonparametric.lowess(clip_y, clip_x, frac=frac, it=0,
        #                             delta=0.0, is_sorted=True)
        le = lowess(clip_y, clip_x, frac=frac, it=0, delta=0.0, is_sorted=True)
        """
        while (stop == False):
            frac = smooth*u.nm/range_x 
            le = sm.nonparametric.lowess(clip_y, clip_x, frac=frac, it=0,
                                         delta=0.0, is_sorted=True)
            cont_y = np.interp(clip_x, le[:, 0], le[:, 1])
            norm_y = clip_y / cont_y
            clip_y = sigmaclip(norm_y, low=3.0, high=3.0)[0]
            #clip_y = sigmaclip(norm_y, low=3.0, high=100.0)[0]            
            clip_x = clip_x[np.isin(norm_y, clip_y)]
            cont_y = cont_y[np.isin(norm_y, clip_y)]
            clip_y = clip_y * cont_y
            stop = (len(clip_y) == len(norm_y)) #or i > 1
        """

        self._clip_x = clip_x
        self._clip_y = clip_y
        self._cont.y = np.interp(self._cont.x, le[:, 0], le[:, 1])
        #self._cont.y = np.interp(self._cont.x, clip_x, cont_y)
        #self._cont.y = spl(self._cont.x)
        
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
        ### WORKS ONLY ON LOGARITHMICALLY SPACED SPECTRA!
        
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



    def convolve(self, col='y', prof=None, gauss_sigma=20, convert=True,
                 mode='same'):
        """Convolve a spectrum with a profile using FFT transform
        
        The profile must have the same length of the column X of the spectrum.
        If no profile @a prof is provided, a normalized Gaussian profile with 
        @a gauss_sigma is applied."""

        conv = copy.deepcopy(self)
        conv_col = getattr(conv, col)
        if convert == True:
            conv.todo_convert_logx(xunit=u.km/u.s)        
        if prof is None:
            par = np.stack([[np.median(conv.x.value)], [gauss_sigma], [1]])
            par = np.ndarray.flatten(par, order='F')
            prof = many_gauss(conv.x.value, *par)
            if (len(prof) % 2 == 0):
                prof = prof[:-1]
            prof = prof / np.sum(prof)
        conv.y = fftconvolve(conv_col, prof, mode=mode)
        if convert == True:
            conv.todo_convert_logx(xunit=self.x.unit)
        return conv
        
    def deredden(self, A_v, model='od94'):
        extFactor = extinction.reddening(self._t['X'], A_v, model=model)
        self._t['Y']  *= extFactor
        self._t['DY'] *= extFactor
    
    def extract(self, xmin=None, xmax=None, prox=False, forest=[], zem=[],
                prox_vel=[], line=None):
        """ Extract a region of a spectrum """

        self._orig = dc(self)
        #if (prox != (zem == [])):
        #    raise Exception("Forest name and emission redshift must be "
        #                    "provided together.")
        if ((forest == []) != (zem == [])):
            raise Exception("Forest name and emission redshift must be "
                            "provided together.")
        #if ((forest != []) == prox):
        #    raise Exception("Please choose either forest or proximity.")

        reg = dc(self)
        if (prox == False):

            # Check this part
            if ((forest != []) and (prox_vel == [])):  
                prox_vel = 7e3 * (u.km/u.s)
            if (forest == 'Ly'):
                xmin = 0*u.nm#dict_wave['Ly_lim'] * (1 + zem + prox_vel / c)
                xmax = dict_wave['Ly_a'] * (1 + zem) * (1 - prox_vel / c)
            if (forest == 'CIV'):
                xmin = dict_wave['Ly_a'] * (1 + zem) * (1 + prox_vel / c)
                xmax = dict_wave['CIV_1550'] * (1 + zem) * (1 - prox_vel / c)
            if (forest == 'MgII'):
                xmin = dict_wave['Ly_a'] * (1 + zem) * (1 + prox_vel / c)
                xmax = dict_wave['MgII_2803'] * (1 + zem) * (1 - prox_vel / c)
                

            lt_xmin = []
            gt_xmax = []        
            if (xmin is not None):
                lt_xmin = np.where(reg.x < xmin)[0]
            if (xmax is not None):
                gt_xmax = np.where(reg.x > xmax)[0]
                
            where = np.hstack((lt_xmin, gt_xmax))

        else:

            ion = dict_doubl['Ly']
            wave_list = np.array([dict_wave[ion[i]].value \
                                  for i in range(len(ion))])
            xmin = wave_list * (1 + zem) * (1 - prox_vel/c)
            xmax = wave_list * (1 + zem)
            where = np.full(len(reg.x), True) 
            for i in range(len(xmin)):
                w = np.logical_and(reg.x.value > xmin[i], reg.x.value < xmax[i])
                where[w] = False
        reg._t.remove_rows(where)

        # Check this part
        if (hasattr(reg, '_cont')):
            reg._cont._t.remove_rows(where)
                
            
        return reg

    def find_extrema(self):
        """Find the extrema in a spectrum and save them as a spectrum"""

        #min_idx = []
        #max_idx = []
        #extr_idx = []
        if (len(self.t) > 0):
            min_idx = np.hstack(argrelmin(self.y))
            max_idx = np.hstack(argrelmax(self.y))
            extr_idx = np.sort(np.append(min_idx, max_idx))
            minima = self.from_table(self._t[min_idx])
            maxima = self.from_table(self._t[max_idx])
            extr = self.from_table(self._t[extr_idx])
        else:
            minima = self.from_table(self._t)
            maxima = self.from_table(self._t)
            extr = self.from_table(self._t)
        return minima, maxima, extr
        
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

        x = []
        y = []
        xmin = []
        xmax = []
        dy = None
        good = None
        resol = None
        xunit = None
        yunit = None
        meta = None
        if (len(table) > 0):
            x = table['X']
            y = table['Y']            
            xmin = table['XMIN']
            xmax = table['XMAX']
            dy = table['DY']
            group = table['GROUP']
            resol = table['RESOL']
            xunit = x.unit
            yunit = y.unit
            dx = 0.5 * (xmax - xmin)
            
            c1 = np.argwhere(y > 0)
            c2 = np.argwhere(dy > 0)
            igood = np.intersect1d(c1, c2)

            good = np.repeat(-1, len(x))
            good[igood] = 1

        spec = Spec1D(x, y, dy=dy, xmin=xmin, xmax=xmax, xunit=xunit, 
                      yunit=yunit, group=good, resol=resol, meta=meta)
        return spec
        
    def mask_col(self, col):
        """ Mask columns """
        
        if self._use_good:
            ret = self._t[col].quantity[self._igood]
        else:
            ret = self._t[col].quantity
        null = np.argwhere(self._t[col].mask)
        if null.size > 0:
            ret[null] = float('nan')
        return ret

    def model(self, prof=True, psf=True):
        model = self._norm_guess[0]
        param = self._norm_guess[1]
        
        if (prof == True):
            model = model * self._prof_guess[0]
            param.update(self._prof_guess[1])
            
        if (psf == True):
            model = lmc(model, self._psf[0], convolve)
            param.update(self._psf[1])

        self._guess = (model, param)
        return self._guess

    def plot(self, ax=None):
        spec = self
        if ax == None:
            fig = plt.figure(figsize=(10,4))
            fig.canvas.set_window_title("Spectrum")
            grid = gs(1, 1)
            ax = fig.add_subplot(grid[:, :])
        #ax_z = ax.twiny()
        ax.set_xlabel("Wavelength [" + str(spec.x.unit) + "]")
        ax.set_ylabel("Flux [" + str(spec.y.unit) + "]")
        #ax_z.set_xlabel("Lyman-alpha redshift")
        ax.plot(spec.x, spec.y, lw=1.0)
        #ax_z.plot(spec.x/dict_wave['Ly_a']-1, spec.y, lw=1.0)
        #ax.plot(spec.x, spec.dy, c='r', lw=1.0)
        if (hasattr(self, '_cont')):
            where = np.where(self.y != self._cont.y)
            #print(self.y, self._cont.y)
            #print(where)
            #ax.plot(self._cont.x[where], self._cont.y[where], c='y')
            ax.plot(self.x, self._cont.y, c='y')

        #if block is False:
        #    plt.ion()
        #    plt.draw()
        #fig.tight_layout()
        if ax == None:
            grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.9])
            grid.update(wspace=0.2, hspace=0.0)
            plt.show()

    def prof_auto2(self, model, ion, logN_range=[12.0], b_range=[8] * u.km/u.s,
                   width=0.03 * u.nm, plot=False, verbose=False):

        loop = [(logN, b) for logN in logN_range for b in b_range]
        match = np.array([])
        minima_x = np.array([])
        minima_redchi = np.array([])
        for (logN, b) in loop:
            N = np.power(10, logN) / u.cm**2
            prof = model.prof_mult(self, ion, b=b, N=N)
            redchi = np.zeros((len(self.x), len(prof)))
            # Not working with Python 2.7
            #if (verbose == True):
                #print("LogN = %.2f, b = %.2f km/s: scanning" % (logN, b.value),
                #      end=" ", flush=True)
            for p in range(len(prof)):
                if (p == 0):
                    redchi[:, p] = self.prof_scan_new(prof[p], ion,
                                                      verbose=verbose)
                else:
                    redchi[:, p] = self.prof_scan_new(prof[p], ion)

            where = np.full(len(self.x), False)
            for x in range(len(self.x)):
                if (np.array_equal(np.sort(redchi[x, :]), redchi[x, :]) \
                    and (redchi[x, 0] < redchi[x, 1] * 0.95) \
                    and (redchi[x, 0] < 10.0)):                            
                    where[x] = True

                    
            minima = None
            redchi_sel = Spec1D(self.x[where], redchi[:, 0][where],
                                xunit=self.x.unit)
            minima, maxima, extr = redchi_sel.find_extrema()
            minima_x = np.append(minima_x, minima.x)
            minima_redchi = np.append(minima_redchi,
                                      np.interp(minima.x, self.x, redchi[:, 0]))

            if (verbose == True):
                print("%i new matches found, %i in total." \
                      % (len(minima.x), len(np.unique(minima_x))))
                
            if (plot == True):
                self.plot(block=False)
                
                fig = plt.figure(figsize=(10,4))
                fig.canvas.set_window_title("Spectrum")
                grid = gs(1, 1)
                ax = fig.add_subplot(grid[:, :])
                grid.tight_layout(fig, rect=[0.02, 0.02, 1, 0.97])
                ax.set_xlabel("Wavelength [" + str(self.x.unit) + "]")
                ax.set_ylabel("Reduced chi squared")
                ax.set_yscale("log", nonposy='clip')
                ax.plot(self.x, redchi[:, 0], c='black', lw=1.0)
                #ax.plot(self.x[where], redchi[:, 0][where], c='black', lw=2.0)
                ax.plot(self.x, redchi[:, 1], c='r', lw=1.0)
                ax.plot(self.x, redchi[:, 2], c='g', lw=1.0)
                ax.plot(self.x, redchi[:, 3], c='b', lw=1.0)
                if (len(minima.t) > 0):
                    ax.scatter(minima.x, minima.y, c='b')
                plt.show()    

        x = np.array([])
        for p in range(len(ion)):
            x = np.append(x, np.asarray(minima.x) \
                          * dict_wave[ion[p]] / dict_wave[ion[0]]) 

        # Improve...
        xmin = x - width.value
        xmax = x + width.value
            
        y = np.interp(np.asarray(x), self.x.value, self.y.value) \
                 * self.y.unit
        dy = np.interp(np.asarray(x), self.x.value, self.dy.value) \
                  * self.y.unit
        from . import Line
        line = Line(xmin=xmin, xmax=xmax, x=x, y=y, dy=dy, xunit=self.x.unit,
                    spec=self)
        return line

    def prof_auto(self, model, ion, logN_range=[12.0], b_range=[8] * u.km/u.s,
                  width=0.03 * u.nm, plot=False, verbose=False):

        loop = [(logN, b) for logN in logN_range for b in b_range]
        match = np.array([])
        minima_x = np.array([])
        minima_redchi = np.array([])
        for (logN, b) in loop:
            N = np.power(10, logN) / u.cm**2

            wave_list = np.array([dict_wave[ion[i]].value \
                                  for i in range(len(ion))]) 
            ion = np.array(ion)[np.argsort(wave_list)]
            prof = model.prof_mult(self, ion, b=b, N=N, plot=plot)
            #print(prof.x)
            redchi = np.zeros((len(self.x), len(prof)))
            # Not working with Python 2.7
            #if (verbose == True):
            #    print("LogN = %.2f, b = %.2f km/s: scanning" % (logN, b.value),
            #          end=" ", flush=True)

            for p in range(len(prof)):
                if (p == 0):
                    redchi[:, p] = self.prof_scan_new(prof[p], ion,
                                                      verbose=verbose)
                else:
                    redchi[:, p] = self.prof_scan_new(prof[p], ion)

            where = np.full(len(self.x), False)
            for x in range(len(self.x)):
                #if ((redchi[x, 0] == np.min(redchi[x, :])) \
                #    and (redchi[x, 1] == np.max(redchi[x, :])) \
                #    and (redchi[x, 2] < redchi[x, 3]) \
                #    and (redchi[x, 0] < 10.0)):
                if (np.array_equal(np.sort(redchi[x, :]), redchi[x, :]) \
                    #and (redchi[x, 0] < redchi[x, 1] * 0.95) \
                    and (redchi[x, 0] < redchi[x, 1] * 0.95) \
                    #and (redchi[x, -2] < redchi[x, -1] * 0.95) \
                    #and (redchi[x, 0] < 10.0) \
                    #and (redchi[x, 0] < redchi[x, 1]) \
                    ):
                    where[x] = True

            minima = None
            #if (np.sum(where) > 0):
            redchi_sel = Spec1D(self.x[where], redchi[:, 0][where],
                                xunit=self.x.unit)
            minima, maxima, extr = redchi_sel.find_extrema()
            minima_x = np.append(minima_x, minima.x)
            minima_redchi = np.append(minima_redchi,
                                      np.interp(minima.x, self.x, redchi[:, 0]))

            #print(minima_x, minima_redchi)
            
            if (verbose == True):
                print("%i new matches found, %i in total." \
                      % (len(minima.x), len(np.unique(minima_x))))
                
            if (plot == True):
                self.plot(block=False)
                
                fig = plt.figure(figsize=(10,4))
                fig.canvas.set_window_title("Spectrum")
                grid = gs(1, 1)
                ax = fig.add_subplot(grid[:, :])
                grid.tight_layout(fig, rect=[0.02, 0.02, 1, 0.97])
                ax.set_xlabel("Wavelength [" + str(self.x.unit) + "]")
                ax.set_ylabel("Reduced chi squared")
                ax.set_yscale("log", nonposy='clip')
                ax.plot(self.x, redchi[:, 0], c='black', lw=1.0)
                #ax.plot(self.x[where], redchi[:, 0][where], c='black', lw=2.0)
                ax.plot(self.x, redchi[:, 1], c='r', lw=1.0)
                ax.plot(self.x, redchi[:, 2], c='g', lw=1.0)
                ax.plot(self.x, redchi[:, 3], c='b', lw=1.0)
                if (len(minima.t) > 0):
                    ax.scatter(minima.x, minima.y, c='b')
                plt.show()    

            #print(np.asarray(minima.x))

        #print(mimima_x, minima_redchi)
        x = np.array([])
        for p in range(len(ion)):
            x = np.append(x, np.asarray(minima.x) \
                          * dict_wave[ion[p]] / dict_wave[ion[0]]) 

        # Improve...
        xmin = x - width.value
        xmax = x + width.value
            
        y = np.interp(np.asarray(x), self.x.value, self.y.value) \
                 * self.y.unit
        dy = np.interp(np.asarray(x), self.x.value, self.dy.value) \
                  * self.y.unit
        from . import Line
        line = Line(xmin=xmin, xmax=xmax, x=x, y=y, dy=dy, xunit=self.x.unit,
                    spec=self)
        return line

    #return minima
        #cond = np.logical_and(
        #for logN in logN_range
    
    def prof_merge(self, prof, ion=['Ly_a'], z=None):
        """ Merge a line profile into a spectrum at a random position """

        if (z == None):
            temp = dc(self)
            temp.to_z([ion[-1]])
            zmin = np.min(temp.x)
            temp = dc(self)
            temp.to_z([ion[0]])
            zmax = np.max(temp.x)
            z = random.uniform(zmin, zmax)
        
        x_prof = prof.x * (1 + z)
        xmin = np.min(x_prof)
        xmax = np.max(x_prof)
        chunk = np.logical_and((self.x.value > xmin.value),
                               (self.x.value < xmax.value))
        y = np.interp(self.x[chunk], x_prof, prof.y) * self.y.unit
        self.y[chunk] = self.y[chunk] * y / self._cont.y[chunk]
        
        return z
    
    def prof_scan(self, ion, line=None, width=0.02 * u.nm, verbose=True):
        """ Scan a spectrum with a line profile to find matches """

        
        
        i_last = 0
        z_temp = np.array([])
        z_match = np.array([])

        spec_temp = dc(self)
        spec_temp.to_z([ion[0]])

        xstart = np.min(self.x).value
        xend = np.max(self.x).value
        null = np.ones(len(self.x))

        redchi_thr = 0.7
        start = time.time()

        self._redchi = np.zeros(len(self.x))
        for i in range(len(self.x)):
            # Not working with Python 2.7
            #if ((i == 999) and (verbose==True)): 
            #    print("Scanning (foreseen running time: %.0f s)..." \
            #          % ((time.time() - start) * len(self.x) * 1e-3), end=" ",
            #          flush=True)
            for prof_i in prof:
                mod_x = prof_i.x * (1 + spec_temp.x[i])
                xmin = np.min(mod_x).value
                xmax = np.max(mod_x).value
                if (xmax < xend):
                    #chunk[:, i] = np.logical_and((self.x.value > xmin),
                    #                             (self.x.value < xmax))
                    chunk_i = np.logical_and((self.x.value > xmin),
                                                 (self.x.value < xmax))
                else:
                    chunk_i = np.full(len(self.x), False)
                    
                #chunk_i = chunk[:, i]
                x = self.x[chunk_i].value
                y = self.y[chunk_i].value
                dy = self.dy[chunk_i].value
                mod_y = np.interp(x, mod_x, prof_i.y) * prof_i.y.unit
                null_y = self._cont.y[chunk_i].value
                if (len(x) > 0):
                    redchi_mod = redchi(x, y, dy, mod_y)

                    #"""
                    redchi_0 = redchi(x, y, dy, null_y)

                    #Check this
                    #if ((redchi_mod < redchi_thr * redchi_0) \
                        #    and (redchi_mod < 2.0)):
                    #if (redchi_mod < redchi_thr):
                    if (redchi_mod < redchi_thr * redchi_0):
                        z_temp = np.append(z_temp, spec_temp.x[i])
                        i_last = i
                    else:
                        if (i == i_last+1):
                            z_match = np.append(z_match, np.mean(z_temp))
                            z_temp = np.array([])
                    #"""
                self._redchi[i] = redchi_mod

        z_match = z_match[~np.isnan(z_match)]
        if (verbose == True):
            print("%i matches found." % len(z_match))
        
        x = np.array((1.0 + np.asarray(z_match)) * dict_wave[ion[0]])
        x = np.append(x, (1.0 + np.asarray(z_match)) * dict_wave[ion[1]]) \
            #* self.x.unit

        # Improve...
        xmin = x - width.value
        xmax = x + width.value
            
        y = np.interp(np.asarray(x), self.x.value, self.y.value) \
                 * self.y.unit
        dy = np.interp(np.asarray(x), self.x.value, self.dy.value) \
                  * self.y.unit
        from . import Line
        line = Line(xmin=xmin, xmax=xmax, x=x, y=y, dy=dy, xunit=self.x.unit,
                    spec=self)
        return line
    
    def prof_scan_new(self, prof, ion, width=0.02 * u.nm, verbose=False):
        """ Scan a spectrum with a line profile to find matches """
        
        i_last = 0
        z_temp = np.array([])
        z_match = np.array([])

        spec_temp = dc(self)
        spec_temp.to_z([ion[0]])

        xstart = np.min(self.x).value
        xend = np.max(self.x).value
        null = np.ones(len(self.x))

        start = time.time()

        #for i in range(len(ion)):
            
        
        self._redchi = np.zeros(len(self.x))
        
        for i in range(len(self.x)):
        #for i in range(999):
            # Not working with Python 2.7
            #if ((i == 999) and (verbose==True)): 
            #    print("(foreseen running time: %.0f s)..." \
            #          % ((time.time() - start) * len(self.x) * (len(ion) + 2) \
            #             * 1e-3), end=" ", flush=True)
            mod_x = prof.x * (1 + spec_temp.x[i])


            # To determine the regions where the model is defined, look for
            # jumps in the x array
            mod_xmin = prof.xmin * (1 + spec_temp.x[i])
            mod_xmax = prof.xmax * (1 + spec_temp.x[i])
            #print(prof.xmin, prof.xmax)
            mod_xmin_l = np.append(mod_xmin, float("inf"))
            mod_xmax_l = np.append(0, mod_xmax)
            #print(mod_xmin, mod_xmax)
            where = np.abs(mod_xmin_l.value-mod_xmax_l.value) \
                    > (mod_xmax[1].value-mod_xmin[0].value)
            xmin = mod_xmin_l[where][::2].value
            xmax = mod_xmax_l[where][1::2].value
            #print(mod_xmin_l[where], mod_xmax_l[where])
            #print(xmin, xmax)
            
            #print(prof.t)
            #xmin = np.min(mod_x).value
            #xmax = np.max(mod_x).value
            #print(np.min(mod_x).value, np.max(mod_x).value)
            
            #"""
            chunk_i = np.full(len(self.x), False)
            for j in range(np.size(xmin)):
                if (xmax[j] < xend):
                    #chunk[:, i] = np.logical_and((self.x.value > xmin),
                    #                             (self.x.value < xmax))
                    chunk_i = np.logical_or(chunk_i,
                                  np.logical_and((self.x.value > xmin[j]),
                                                 (self.x.value < xmax[j])))
                else:
                    chunk_i = np.full(len(self.x), False)
            #"""        
            #if (xmax < xend):
            #    #chunk[:, i] = np.logical_and((self.x.value > xmin),
            #    #                             (self.x.value < xmax))
            #    chunk_i = np.logical_and((self.x.value > xmin),
            #                                 (self.x.value < xmax))
            #else:
            #    chunk_i = np.full(len(self.x), False)

            #chunk_i = np.logical_and(
            
            #chunk_i = chunk[:, i]
            x = self.x[chunk_i].value
            y = self.y[chunk_i].value/self._cont.y[chunk_i].value
            dy = self.dy[chunk_i].value/self._cont.y[chunk_i].value
            mod_y = np.interp(x, mod_x, prof.y) #* prof.y.unit
            if (len(x) > 0):
                redchi_mod = redchi_f(x, y, dy, mod_y)
            else:
                redchi_mod = 'nan'
            self._redchi[i] = redchi_mod

            #plt.plot(x, y)
            #plt.plot(x, dy)
            #plt.plot(x, mod_y)
            #plt.show()
            
        return self._redchi

    def redchi(self, model_param=None, nvarys=0, chunk=None):
        if (hasattr(self, '_spec')):
            spec = self._spec
        else:
            spec = self
        model = model_param[0]
        param = model_param[1]
        #param.pretty_print()
        if (chunk is None):
            x = spec.x
            y = spec.y      
            dy = spec.dy
        else:
            x = spec.x[chunk]
            y = spec.y[chunk]        
            dy = spec.dy[chunk]
        ndata = len(x)
        mod = model.eval(param, x=x.value)
        ret = np.sum(((mod-y.value)/dy.value)**2) / (ndata-nvarys)
        return ret
    
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
            [fits.Column(name='XMIN', format='E', array=self.xmin,
                         unit=self.xmin.unit.name),
             fits.Column(name='XMAX', format='E', array=self.xmax,
                         unit=self.xmax.unit.name),
             fits.Column(name='X', format='E', array=self.x,
                         unit=self.x.unit.name),
             fits.Column(name='Y', format='E', array=self.y,
                         unit=self.y.unit.to_string()),
             fits.Column(name='DY', format='E', array=self.dy,
                         unit=self.dy.unit.to_string()),
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
        self.x = np.unique(x.value) * x.unit
        self.xmin = np.unique(xmin) * x.unit
        self.xmax = np.unique(xmax) * x.unit

    def to_z(self, ion):
        """ Convert line wavelengths to redshifts using a reference """

        if (len(ion) == 1):
            ion = np.full(len(self.t), ion[0])

        if ((len(ion) > 1) and (len(ion) != len(self.t))):
            raise Exception("Ion list and table must have the same length.")
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

