from . import Spec1D, Model
from .utils import convolve, dict_wave, savitzky_golay, voigt_def
from astropy.table import Column, Table
from astropy.stats import sigma_clip
from astropy import units as u
from astropy.io import fits
from copy import deepcopy as dc
import inspect
from lmfit import CompositeModel as lmc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)
from scipy.interpolate import UnivariateSpline
import scipy.fftpack as fft
from scipy.signal import argrelmin, argrelmax, fftconvolve, savgol_filter, wiener
from scipy.stats import sigmaclip
#from statsmodels.nonparametric.kernel_regression import KernelReg
#import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess
import sys
import time
import warnings

yunit = u.erg / (u.Angstrom * u.cm**2 * u.s)


class Line(Spec1D):

    def __init__(self,
                 spec=None,
                 x=[],
                 y=[],
                 xmin=[],
                 xmax=[],
                 dy=[],
                 group=[],
                 resol=[],
                 xunit=None,
                 yunit=None,
                 meta=None,
                 dtype=float):
        """ Constructor for the Line class """
        
        # Exceptions
        if ((x == []) != (y == [])):
            raise Exception("X and Y must be provided together.")
        if ((xmin == []) != (xmax == [])):
            raise Exception("XMIN and XMAX must be provided together.")
        sumlen = len(x) + len(y) + len(xmin) + len(xmax) + len(dy)
        if ((x != []) and (sumlen % len(x) != 0)):
            raise Exception("Data arrays must have the same length.")
        if ((xmin != []) and (sumlen % len(xmin) != 0)):
            raise Exception("Data arrays must have the same length.")

        # Warnings
        if ((spec is None) and (x is [])):
            warnings.warn("No spectrum or data provided.")

        # Spectrum
        self._spec = None
        if (spec is not None):
            self._spec = dc(spec)
            if (hasattr(spec, '_precont')):
                self._precont = dc(spec._precont)
            if (hasattr(spec, '_cont')):
                self._cont = dc(spec._cont)            
                
        # Line list
        data = ()
        #if (x != []):
        if (np.size(x)>0):
            col_x = Column(np.asarray(dc(x), dtype=dtype), name='X')
            col_y = Column(np.asarray(dc(y), dtype=dtype), name='Y')
            data = (col_x, col_y)
            if xunit is None:
                try:
                    xunit = x.unit
                except:
                    raise Exception("X unit not provided.")    
            if yunit is None:
                try:
                    yunit = y.unit
                except:
                    raise Exception("Y unit not provided.")
        #if (xmin != []):
        if (np.size(xmin) > 0):
            col_xmin = Column(np.asarray(dc(xmin), dtype=dtype), name='XMIN')
            col_xmax = Column(np.asarray(dc(xmax), dtype=dtype), name='XMAX')
            data += (col_xmin, col_xmax)
        #if (dy != []):
        if (np.size(dy) > 0):
            col_dy = Column(np.asarray(dc(dy), dtype=dtype), name='DY')
            data += (col_dy,)

        # Needed for backwards compatibility with Spec1DCont
        if (group == []):
            group = np.ones(len(x))
        col_group = Column(np.asarray(dc(group), dtype=dtype), name='GROUP')
        data += (col_group,)
        if (resol == []):
            resol = np.zeros(len(x))
        col_resol = Column(np.asarray(dc(resol), dtype=dtype), name='RESOL')
        data += (col_resol,)
        
        if data is ():
            data = None
            
        if (meta is None):
            meta = {}

        self._t = Table(data=data, masked=True, meta=meta)
        if (x != []):
            self._t['X'].unit = xunit
            self._t['Y'].unit = yunit
        if (xmin != []):
            self._t['XMIN'].unit = xunit
            self._t['XMAX'].unit = xunit
        if (dy != []):
            self._t['DY'].unit = yunit

        self._ion = np.full(len(self._t), 'Ly_a')           


        self._use_good = False

        
# Properties
        
    @property
    def dy(self):
        return self.mask_col('DY')

    @dy.setter
    def dy(self, value):
        if self._use_good:
            self._t['DY'][self._igood] = np.asarray(value, dtype=float)
        else:
            self._t['DY'] = np.asarray(value, dtype=float)
        try:
            self._t['DY'].unit = value.unit
        except:
            raise Exception("DY unit not provided.")

    @property
    def t(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t
        
    @property
    def spec(self):
        return self._spec

    @spec.setter
    def spec(self, value):
        if isinstance(value, Spec1D):
            self._spec = value
        else:
            raise Exception("Spectrum has a wrong format.")
    
    @property
    def use_good(self):
        return self._use_good

    @use_good.setter
    def use_good(self, value):
        self._use_good = value
        if self._use_good:
            self._igood = np.argwhere(self._t['GROUP'] >= 0)

    @property
    def x(self):
        return self.mask_col('X')

    @x.setter
    def x(self, value):
        if self._use_good:
            self._t['X'][self._igood] = np.asarray(value, dtype=float)
        else:
            self._t['X'] = np.asarray(value, dtype=float)
        try:
            self._t['X'].unit = value.unit
        except:
            raise Exception("X unit not provided.")

    @property
    def xmin(self):
        return self.mask_col('XMIN')

    @xmin.setter
    def xmin(self, value):
        if self._use_good:
            self._t['XMIN'][self._igood] = np.asarray(value, dtype=float)
        else:
            self._t['XMIN'] = np.asarray(value, dtype=float)
        try:
            self._t['XMIN'].unit = value.unit
        except:
            raise Exception("XMIN unit not provided.")

    @property
    def xmax(self):
        return self.mask_col('XMAX')

    @xmax.setter
    def xmax(self, value):
        if self._use_good:
            self._t['XMAX'][self._igood] = np.asarray(value, dtype=float)
        else:
            self._t['XMAX'] = np.asarray(value, dtype=float)
        try:
            self._t['XMAX'].unit = value.unit
        except:
            raise Exception("XMAX unit not provided.")

    @property
    def y(self):
        return self.mask_col('Y')

    @y.setter
    def y(self, value):
        if self._use_good:
            self._t['Y'][self._igood] = np.asarray(value, dtype=float)
        else:
            self._t['Y'] = np.asarray(value, dtype=float)
        try:
            self._t['Y'].unit = value.unit
        except:
            raise Exception("Y unit not provided.")

        
# Methods
    
    """
    def chunk(self, x=None, line=None, single=False):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x.value-x.value) \
                            == abs(self.x.value-x.value).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        iter = range(len(self._t))
        if (hasattr(self._spec, '_orig')):
            spec = dc(self._spec._orig)
        else:
            spec = dc(self._spec)
        sel = spec.t['X'] < 0.0
        for row in self.t[self.group(line=line, single=single)[1]]:
            sel = np.logical_or(sel, np.logical_and(
                spec.t['X'] >= row['XMIN'],
                spec.t['X'] <= row['XMAX']))
        if (np.sum(sel) % 2 == 0):
            sel[np.argmax(sel)] = 0
        ret = (line, sel)
        self._chunk = ret
        for c in range(1, len(ret)):
            if (c == 1):
                self._chunk_sum = dc(ret[c])
            else:
                self._chunk_sum += ret[c]
        return ret
    """

    def cont(self, smooth=3.0, flux_corr=1.0, kappa_low=2.0, kappa_high=2.0):
        """ Determine the emission continuum by removing absorption lines """
        
        self._precont = dc(self._spec)
        range_x = np.max(self._spec.x) - np.min(self._spec.x)
        x = self._maxs['X']
        y = self._maxs['Y'] * flux_corr
        self._cont = dc(self._precont)        

        clip_x = x
        clip_y = y
        stop = False
        i = 0


        while (stop == False):
            #plt.plot(self._spec.x, self._spec.y, color='black')
            #plt.scatter(clip_x, clip_y, color='red')
            #plt.show()
            #frac = min(1, smooth * 10/len(x)) 
            #le = sm.nonparametric.lowess(clip_y, clip_x, frac=frac)
            frac = smooth*u.nm/range_x #* len(x)/len(clip_x)
            #le = sm.nonparametric.lowess(clip_y, clip_x, frac=frac, it=0,
            #                             delta=0.0, is_sorted=True)
            le = lowess(clip_y, clip_x, frac=frac, it=0, delta=0.0,
                        is_sorted=True)
            cont_y = np.interp(clip_x, le[:, 0], le[:, 1])
            norm_y = clip_y / cont_y
            #clip_y = sigmaclip(norm_y, low=2.0, high=100.0)[0]
            clip_y = sigmaclip(norm_y, low=kappa_low, high=kappa_high)[0]
            #clip_y = sigmaclip(norm_y, low=2.0, high=100.0)[0]            
            clip_x = clip_x[np.isin(norm_y, clip_y)]
            cont_y = cont_y[np.isin(norm_y, clip_y)]
            clip_y = clip_y * cont_y
            stop = ((len(clip_y) == len(norm_y)) \
                    or (len(clip_y) < 100))


        self._clip_x = clip_x
        self._clip_y = clip_y
        self._cont.t['Y'] = np.interp(self._cont.x, le[:, 0], le[:, 1]) \
                            * self._cont.y.unit

        
    """    
    def corr_resid(self, group, cont_corr):

        xmin = np.min(self.xmin[group[1]])
        xmax = np.max(self.xmax[group[1]])
        fit_where = np.where(np.logical_and(self._resid_fit.x >= xmin,
                                            self._resid_fit.x <= xmax))
        cont_where = np.where(np.logical_and(self._resid_cont.x >= xmin,
                                             self._resid_cont.x <= xmax))
        fit_norm = self._resid_fit.y.value * 0.0
        fit_norm[fit_where] = self._resid_fit.y[fit_where] /\
                              self._resid_fit.dy[fit_where]
        cont_norm = self._resid_cont.y.value * 0.0
        cont_norm[cont_where] = self._resid_cont.y[cont_where] /\
                                self._resid_cont.dy[cont_where]
        cont_pos = self._resid_cont.y.value * 0.0
        cont_pos[cont_where] = self._resid_cont.y[cont_where]
        cont_pos = cont_pos[cont_pos > 0]
        cont_rms = 0.0
        if (len(cont_pos) > 0):
            cont_rms = np.sqrt(np.mean(cont_pos**2))
        
        x_lo = self._resid_fit.x[np.argmin(fit_norm)]
        x_hi = self._resid_cont.x[np.argmax(cont_norm)]
        y_lo = np.min(fit_norm)
        y_hi = np.max(cont_norm)
        if (np.absolute(y_lo) > np.absolute(y_hi) or 1==1):
            y = np.interp(x_lo.value, self._spec.x, self._spec.y) 
            dy = np.interp(x_lo.value, self._spec.x, self._spec.dy)
            self._t.add_row([x_lo, y, xmin, xmax, dy])
            cont_corr = 1.0
        if (cont_rms > 3*np.mean(self._resid_cont.dy.value)):
        #else:
            cont_corr = 1 + cont_rms / np.mean(self._cont.y.value)
        
        # This is actually not needed for un-identified lines (no need to 
        # differentiate between components of the current system and spurious
        # components) but is now used to allow compatibility with the 
        # syst.corr_resid method. Fix that.
        self._noneb = dc(self)
        self._neb = dc(self)
        
        return cont_corr
    """

    def corr_tau(self, N_thres=None):

        tau_norm = 0.0028
        tau_index = 3.45
        logN_1 = 12.0
        logN_2 = 14.0
        logN_3 = 18.0
        logN_4 = 19.0
        logN_index = -1.65
        logN_thres = np.log(N_thres.value)
        fact_0 = 1.0
        fact_1 = 1e14 / np.log(1e14)
        fact_2 = np.log(1e18) / np.log(1e14) * 1e5
        if (logN_thres <= logN_1):
            num = 0
        elif ((logN_thres > logN_1) and (logN_thres <= logN_2)):
            num = (pow(10, logN_thres * (2.0 + logN_index)) \
                   - pow(10, logN_1 * (2.0 + logN_index))) / (2.0 + logN_index)
        elif ((logN_thres > logN_2) and (logN_thres <= logN_3)): 
            num = (pow(10, logN_2 * (2.0 + logN_index)) \
                   - pow(10, logN_1 * (2.0 + logN_index))) / (2.0 + logN_index) 
            num = num + fact_1 \
                  * (pow(10, logN_thres * (1.0 + logN_index))) \
                  /(1.0 + logN_index) \
                  * (np.log(pow(10, logN_thres)) - 1.0 / (1.0 + logN_index))
            num = num - fact_1 \
                  * (pow(10, logN_2 * (1.0 + logN_index))) \
                  / (1.0 + logN_index) \
                  * (np.log(pow(10, logN_2)) - 1.0 / (1.0 + logN_index))
        elif ((logN_thres > logN_3) and (logN_thres <= logN_4)): 
            num = (pow(10, logN_2 * (2.0 + logN_index)) \
                   - pow(10, logN_1 * (2.0 + logN_index))) / (2.0 + logN_index)
            num = num + fact_1 \
                  * (pow(10, logN_3 * (1.0 + logN_index))) \
                  / (1.0 + logN_index) \
                  * (np.log(pow(10, logN_3)) - 1.0 / (1.0 + logN_index))
            num = num - fact_1 \
                  * (pow(10, logN_2 * (1.0 + logN_index))) \
                  / (1.0 + logN_index) \
                  * (np.log(pow(10, logN_2)) - 1.0 / (1.0 + logN_index))
            num = num + fact_2 \
                  * (pow(10, logN_thres * (1.5 + logN_index)) \
                  - pow(10, logN_3   * (1.5 + logN_index))) / (1.5 + logN_index)
        
        if (logN_thres > logN_4): 
            num = 1
        
        den = (pow(10, logN_2 * (2.0 + logN_index)) \
               - pow(10, logN_1 * (2.0 + logN_index))) / (2.0 + logN_index)
        den = den + fact_1 \
              * (pow(10, logN_3 * (1.0 + logN_index))) / (1.0 + logN_index) \
              * (np.log(pow(10, logN_3)) - 1.0 / (1.0 + logN_index))
        den = den - fact_1 \
              * (pow(10, logN_2 * (1.0 + logN_index))) / (1.0 + logN_index) \
              * (np.log(pow(10, logN_2)) - 1.0 / (1.0 + logN_index))

        z = self._spec.x / dict_wave['Ly_a'] - 1 
        flux_corr = np.exp(tau_norm * pow(1+z, tau_index) * num/den)
        return flux_corr
        
    def find_special(self, mode='abs', diff='max', kappa=3.0, 
                     sigma_min=5.0, sigma_max=100.0):
    
        x = np.array([])
        y = np.array([])
        xmin = np.array([])
        xmax = np.array([])
        dy = np.array([])
        #for log_sigma in np.arange(0.6, 2.0, 0.4):
        #for log_sigma in np.arange(np.log10(sigma_max), np.log10(sigma_min), -0.1):
        for log_sigma in np.arange(np.log10(sigma_min), np.log10(sigma_max), 0.5):
            sigma = 10**log_sigma
        #for sigma in np.arange(sigma_min, sigma_max, 10.0):
            (x_t, y_t, xmin_t, xmax_t, dy_t) = \
                self.find(mode, diff, True, kappa, sigma)
            w = np.zeros(len(x_t), dtype=bool)
            #for min, max in zip(xmin, xmax):
            #    w += np.logical_and(x_t>min, x_t<max)
            #print len(w)
            if (len(w) > 0):
                w_not = np.logical_not(w)
                x = np.append(x, x_t[w_not])
                y = np.append(y, y_t[w_not])
                xmin = np.append(xmin, xmin_t[w_not])
                xmax = np.append(xmax, xmax_t[w_not])
                dy = np.append(dy, dy_t[w_not])
            else:
                x = np.append(x, 0)
                y = np.append(y, 0)
                xmin = np.append(xmin, 0)
                xmax = np.append(xmax, 0)
                dy = np.append(dy, 0)
            x, ind = np.unique(x, return_index=True)
            y = y[ind]
            xmin = xmin[ind]
            xmax = xmax[ind]
            dy = dy[ind]

        line = Line(self.spec, x=x, y=y, xmin=xmin, xmax=xmax, dy=dy, 
                    xunit=u.nm, yunit=yunit)
        self.__dict__.update(line.__dict__)
        self._t.sort('X')   

    def find(self, mode='abs', diff='max', convert=True, kappa=3.0, sigma=10.0):
        """ Find lines """

        if (self._spec is None):
            raise Exception("Spectrum not provided.")

        spec = dc(self._spec)

        # Convolve the 1D spectrum with a gaussian filter
        conv = spec.convolve(gauss_sigma=sigma, convert=convert)
        
        #conv = dc(spec)
        #conv._t['Y'] = savgol_filter(spec._t['Y'], int(sigma)*2+1, 4,  mode='nearest')
        #conv._t['Y'] = wiener(spec._t['Y'], int(sigma)*2+1)
        
        self._conv = conv

        
        
        # Find extrema
        mins, maxs, exts = conv.find_extrema()
        self._mins = mins
        self._maxs = maxs
        self._exts = exts

        spec_exts_y = np.interp(exts['X'], spec.t['X'], spec.t['Y']) \
                      * spec.t['Y'].unit
        spec_exts_dy = np.interp(exts['X'], spec.t['X'], spec.t['DY']) \
                       * spec.t['DY'].unit
        
        # Compute the difference between each extremum and its neighbours
        # N.B. Negative fluxes give wrong results! To be fixed 
        #diff_y_left = (exts['Y'][:-2] - exts['Y'][1:-1]) 
        #diff_y_right = (exts['Y'][2:] - exts['Y'][1:-1]) 
        diff_y_left = (spec_exts_y[:-2] - spec_exts_y[1:-1]) 
        diff_y_right = (spec_exts_y[2:] - spec_exts_y[1:-1]) 
        if mode is 'em':
            diff_y_left = -diff_y_left
            diff_y_right = -diff_y_right
        
        # Check if the difference is above threshold
        diff_y_min = np.minimum(diff_y_left, diff_y_right)
        diff_y_max = np.maximum(diff_y_left, diff_y_right)
        thres = exts['DY'][1:-1] * kappa
        if diff is 'min':
            pos = np.greater(diff_y_min, thres)
        else:
            pos = np.greater(diff_y_max, thres)


        """
        x = exts['X'][1:-1][pos]
        #y = exts['Y'][1:-1][pos]
        #dy = exts['DY'][1:-1][pos]
        #y = np.interp(exts['X'][1:-1][pos], spec.t['X'], spec.t['Y']) \
        #              * spec.t['Y'].unit
        #dy = np.interp(exts['X'][1:-1][pos], spec.t['X'], spec.t['DY']) \
        #               * spec.t['DY'].unit
        y = spec_exts_y[1:-1][pos]
        dy = spec_exts_dy[1:-1][pos]
        
        
        xmin = []
        xmax = []        
        if len(x) > 0: 

            #xmin = np.asarray(exts['X'][:-2][pos]) * x.unit
            #xmax = np.asarray(exts['X'][2:][pos]) * x.unit
            
            xmin = (1 - 3*sigma/300000)*x
            xmax = (1 + 3*sigma/300000)*x
        """

        # Select the region of the spectrum between two adjacent maxima
        x = np.array([0])
        xmin = np.array([])
        xmax = np.array([])
        y = np.array([])
        dy = np.array([])
        for lo, up in zip(exts['X'][:-2][pos], exts['X'][2:][pos]):
            sel = np.where(np.logical_and(spec.t['X'] > lo, spec.t['X'] < up))
            cho = np.argmin(spec.t['Y'][sel])
            x = np.append(x, spec.t['X'][sel][cho])
            y = np.append(y, spec.t['Y'][sel][cho])   
            xmin = np.append(xmin, lo)
            xmax = np.append(xmax, up)
            dy = np.append(dy, spec.t['DY'][sel][cho])
        x = x[1:]
        x_unit = spec.t['X'].unit
        y_unit = spec.t['Y'].unit
        x = x * x_unit
        y = y * y_unit
        xmin = xmin * x_unit
        xmax = xmax * x_unit
        dy = dy * y_unit

        self._bound_x = np.setxor1d(xmin, xmax)
        self._bound_y = np.interp(self._bound_x, spec.x, spec.y)
        
        line = Line(self.spec, x=x, y=y, xmin=xmin, xmax=xmax, dy=dy)
        self.__dict__.update(line.__dict__)
        return (x, y, xmin, xmax, dy)
        

    """
    def group(self, x=None, line=None, single=False):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x.value-x.value) \
                            == abs(self.x.value-x.value).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("line is too large.")
        if (single == True):
            sel = np.full(len(self._t), False)
            sel[line] = True
        else:
            iter = range(len(self._t))
            self._t.sort('X')  # This gives a warning on a masked table
            #self._t.group_by('X')
            groups = np.append([0], np.cumsum(self._t['XMAX'][:-1] <
                                                self._t['XMIN'][1:])) 
            sel = np.array([groups[l] == groups[line] for l in iter])
        
        ret = (line, sel)
        self._group = ret
        return ret
    """

    """
    def mask_col(self, col):

        if self._use_good:
            ret = self._t[col].quantity[self._igood]
        else:
            ret = self._t[col].quantity
        null = np.argwhere(self._t[col].mask)
        if null.size > 0:
            ret[null] = float('nan')
        return ret
    """

    """
    def norm(self, group, chunk, value=1.0, vary=False):
        model = Model(self._spec, line=self, group=group, chunk=chunk) 
        norm = model.norm(value, vary)
        if (hasattr(self, '_norm') == False):
            self._norm = dc(self._spec)
        self._norm.y[chunk[1]] = norm[0].eval(
            norm[1], x=self._norm.x[chunk[1]].value) \
            * self._cont.y[chunk[1]].value * self._norm.y[chunk[1]].unit 

        return norm 
    """
    def plot_new(self, ax=None):
        spec = self._spec
        line = self
        if ax == None:
            fig = plt.figure(figsize=(10,4))
            fig.canvas.set_window_title("Lines")
            grid = gs(1, 1)
            ax = fig.add_subplot(grid[:, :])
            ax.plot(spec.x, spec.y, lw=1.0)
        ax.scatter(line.x, line.y, c='r', marker='+', s=100)
        ax.set_xlabel("Wavelength [" + str(line.x.unit) + "]")
        ax.set_ylabel("Flux [" + str(line.y.unit) + "]")
        try:
            ax.plot(spec.x, line._conv.y, linestyle=':')
        except:
            pass

        if ax == None:
            grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.9])
            grid.update(wspace=0.2, hspace=0.0)
            plt.show()

    def plot(self, group=None, chunk=None, figsize=(10,4), block=True,
             **kwargs):
        spec = self._spec
        fig = plt.figure(figsize=figsize)
        fig.canvas.set_window_title("Lines")
        grid = gs(1, 1)
        ax = fig.add_subplot(grid[:, :])
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
        ax.plot(spec.x, spec.y, c='black', lw=1.0)
        ax.plot(spec.x, spec.dy, c='r', lw=1.0)
        ax.set_xlabel("Wavelength [" + str(spec.x.unit) + "]")
        ax.set_ylabel("Flux [" + str(spec.y.unit) + "]")
        ax.scatter(self.x, self.y, c='r', marker='+')
        #if (hasattr(self, '_conv')):
        #    ax.plot(self._conv.x, self._conv.y, c='r')        
        #if (hasattr(self, '_precont')):
        #    ax.plot(self._precont.x, self._precont.y, c='r',
        #            linestyle=':')        
        if (hasattr(self, '_cont')):
            ax.plot(self._cont.x, self._cont.y, c='r', lw=1.0, linestyle=':')
        if (hasattr(self, '_unabs')):
            where = np.where(self._spec.y != self._unabs.y)
            ax.plot(self._unabs.x[where], self._unabs.y[where], c='y', lw=1.0,
                    linestyle=':')
        if (hasattr(self, '_norm')):
            where = np.where(self._spec.y != self._norm.y)
            ax.plot(self._norm.x[where], self._norm.y[where], c='y', lw=1.0,
                    linestyle=':')
        if (chunk is not None):
            ax.plot(spec.x[chunk[1]], spec.y[chunk[1]], c='black', lw=2.0)
            if (hasattr(self, '_voigt')):
                where = np.where(self._spec.y != self._voigt.y)
                ax.plot(self._voigt.x[where], self._voigt.y[where], c='g',
                        lw=1.0, linestyle=':')
        if (hasattr(self, '_cont')):
            where = np.where(self._spec.y != self._cont.y)
            ax.plot(self._cont.x[where], self._cont.y[where], c='y')
        if (hasattr(self, '_fit')):
            where = np.where(self._spec.y != self._fit.y)
            ax.plot(self._fit.x[where], self._fit.y[where], c='g')
        #if (hasattr(self, '_rem')):
        #    where = np.where(self._spec.y != self._rem.y)
        #    ax.plot(self._rem.x[where], self._rem.y[where], c='black',
        #            linestyle=':')
        #if (hasattr(self, '_cont')):
            #ax.scatter(self._extr.x, self._extr.y, c='g')
            #ax.scatter(self._clip_x, self._clip_y, c='r', marker='+')        
        if (group is not None):
            #ax.scatter(self.x[group[0]], self.y[group[0]], c='g', marker='+',
            #           s=100)        
            ax.scatter(self.x[group[1]], self.y[group[1]], c='r', marker='+',
                       s=100)        
        if (hasattr(self, '_fit')):
            fig.suptitle("Reduced chi-squared: %3.2f" % (self._redchi),
                         fontsize=10)
        
        #if block is False:
        #    plt.ion()
        #    plt.draw()
        if block is True:
            plt.show()

    def plot_N(self, N_fit=None, figsize=(4,4), block=True):
        fig = plt.figure(figsize=figsize)
        fig.canvas.set_window_title("Column densities")
        grid = gs(1, 1)
        ax = fig.add_subplot(grid[:, :])
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])

        if (N_fit is None):
            N_fit = self._N_fit
        n, bins, patches = ax.hist(np.log10(N_fit.value), 20)
        if block is True:
            plt.show()
        
    """
    def psf(self, group, chunk, resol):
        
        model = Model(self._spec, line=self, group=group, chunk=chunk)
        psf = model.psf(resol)

        return psf   
    """

    """
    def redchi(self, model_param, nvarys):
        model = model_param[0]
        param = model_param[1]
        x = self._spec.x[self._chunk[1]]
        y = self._spec.y[self._chunk[1]]        
        dy = self._spec.dy[self._chunk[1]]
        ndata = len(x)
        mod = model.eval(param, x=x.value)
        ret = np.sum(((mod-y.value)/dy.value)**2) / (ndata-nvarys)
        return ret
    """

    """
    def unabs(self, group, chunk):
        model = Model(self._spec, line=self, group=group, chunk=chunk) 
        unabs = model.unabs()
        if (hasattr(self, '_unabs') == False):
            self._unabs = dc(self._spec)
        self._unabs.y[chunk[1]] = unabs[0].eval(
            unabs[1], x=self._unabs.x[chunk[1]].value) \
            * self._unabs.y[chunk[1]].unit

        return unabs 
    """

    """
    def voigt(self, group, chunk, z=[], N=[], b=[], btur=[]):

        sumlen = len(z) + len(N) + len(b) + len(btur)
        if ((z != []) and (sumlen % len(z) != 0)):
            raise Exception("Parameter arrays must have the same length.")
        
        model = Model(self._spec, line=self, group=group, chunk=chunk)
        if (z == []):
            z = self.x[group[1]] / dict_wave['Ly_a'] - 1
            if (hasattr(self, '_norm')):
                N = model.N_guess(self._norm)
            else:
                N = model.N_guess(self._unabs)
            if (hasattr(self, '_unabs')):
                cont = self._unabs
            elif (hasattr(self, '_norm')):
                cont = self._norm
            else:
                raise Exception("Continuum not found.")
            N = model.N_guess(cont, ion=self._ion)
            b = np.full(len(self.x[group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(self.x[group[1]]), voigt_def['btur']) \
                   * u.km / u.s
        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)
        voigt = model.voigt(z, N, b, btur, ['Ly_a'])
        self._voigt.y[chunk[1]] = voigt[0].eval(
            voigt[1], x=self._voigt.x[chunk[1]].value) * self._voigt.y.unit
        if (hasattr(self, '_norm')):
            self._voigt.y[chunk[1]] = self._voigt.y[chunk[1]] \
                                      * self._norm.y[chunk[1]].value
        else:
            self._voigt.y[chunk[1]] = self._voigt.y[chunk[1]] \
                                      * self._unabs.y[chunk[1]].value

        self._z_arr = dc(model._z)
        self._N_arr = dc(model._N)
        self._b_arr = dc(model._b)
        self._btur_arr = dc(model._btur)

        return voigt
    """
    def save(self, filename):
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='X', format='E', array=self.x,
                         unit=self.x.unit.name),
             fits.Column(name='XMIN', format='E', array=self.xmin,
                         unit=self.xmin.unit.name),
             fits.Column(name='XMAX', format='E', array=self.xmax,
                         unit=self.xmax.unit.name),
             fits.Column(name='Y', format='E', array=self.y,
                         unit=self.y.unit.to_string()),
             fits.Column(name='DY', format='E', array=self.dy,
                         unit=self.dy.unit.to_string())])
        hdu.writeto(filename, overwrite=True)

