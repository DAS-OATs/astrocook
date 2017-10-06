from . import Spec1D, Model
from .utils import convolve, dict_wave, redchi_thr, savitzky_golay, voigt_def
from astropy.table import Column, Table
from astropy.stats import sigma_clip
from astropy import units as u
from copy import deepcopy as dc
import inspect
from lmfit import CompositeModel as lmc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import scipy.fftpack as fft
from scipy.signal import argrelmin, argrelmax, fftconvolve
from scipy.stats import sigmaclip
from statsmodels.nonparametric.kernel_regression import KernelReg
import statsmodels.api as sm
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
            if (hasattr(spec, '_cont')):
                self._precont = dc(spec._precont)
                self._cont = dc(spec._cont)            
        
        # Line list
        data = ()
        if (x != []):
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
        if (xmin != []):
            col_xmin = Column(np.asarray(dc(xmin), dtype=dtype), name='XMIN')
            col_xmax = Column(np.asarray(dc(xmax), dtype=dtype), name='XMAX')
            data += (col_xmin, col_xmax)
        if (dy != []):
            col_dy = Column(np.asarray(dc(dy), dtype=dtype), name='DY')
            data += (col_dy,)
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
        sel = self._spec.t['X'] < 0.0
        for row in self.t[self.group(line=line, single=single)[1]]:
            sel = np.logical_or(sel, np.logical_and(
                self._spec.t['X'] >= row['XMIN'],
                self._spec.t['X'] <= row['XMAX']))
        if (np.sum(sel) % 2 == 0):
            sel[np.argmax(sel)] = 0
        ret = (line, sel)
        self._chunk = ret
        return ret
    
    def cont(self, smooth=3.0):
        self._precont = dc(self._spec)
        range_x = np.max(self._spec.x) - np.min(self._spec.x)
        x = self._maxima.x
        y = self._maxima.y
        #x = self._extr.x
        #y = self._extr.y
        self._cont = dc(self._precont)        
        clip_x = x
        clip_y = y
        stop = False
        i = 0
        while (stop == False):
            #frac = min(1, smooth * 10/len(x)) 
            #le = sm.nonparametric.lowess(clip_y, clip_x, frac=frac)
            frac = smooth*u.nm/range_x #* len(x)/len(clip_x)
            le = sm.nonparametric.lowess(clip_y, clip_x, frac=frac, it=0,
                                         delta=0.0, is_sorted=True)
            cont_y = np.interp(clip_x, le[:, 0], le[:, 1])
            norm_y = clip_y / cont_y
            #clip_y = sigmaclip(norm_y, low=2.0, high=100.0)[0]
            clip_y = sigmaclip(norm_y, low=3.0, high=3.0)[0]
            clip_y = sigmaclip(norm_y, low=3.0, high=100.0)[0]            
            clip_x = clip_x[np.isin(norm_y, clip_y)]
            cont_y = cont_y[np.isin(norm_y, clip_y)]
            clip_y = clip_y * cont_y
            stop = (len(clip_y) == len(norm_y)) #or i > 1


        self._clip_x = clip_x
        self._clip_y = clip_y
        self._cont.t['Y'] = np.interp(self._cont.x, le[:, 0], le[:, 1]) \
                            * self._cont.y.unit
        #self._cont.y = np.interp(self._cont.x, clip_x, cont_y)
        #self._cont.y = spl(self._cont.x)
    
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
        if (np.absolute(y_lo) > np.absolute(y_hi)):
            y = np.interp(x_lo.value, self._spec.x, self._spec.y) 
            dy = np.interp(x_lo.value, self._spec.x, self._spec.dy)
            self._t.add_row([x_lo, y, xmin, xmax, dy])
            cont_corr = 1.0
        if (cont_rms > 3*np.mean(self._resid_cont.dy.value)):
        #else:
            cont_corr = 1 + cont_rms / np.mean(self._cont.y.value)
        return cont_corr

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
        

    def find(self, mode='abs', diff='max', convert=True, kappa=3.0, sigma=10.0,
             hwidth=2):
        """ Find lines in a spectrum """

        if (self._spec is None):
            raise Exception("Spectrum not provided.")
        
        spec = dc(self._spec)

        # Convolve the 1D spectrum with a gaussian filter
        conv = spec.convolve(gauss_sigma=sigma, convert=convert)
        #conv = dc(spec)
        self._conv = conv
        
        # Find extrema
        minima, maxima, extr = conv.find_extrema()
        self._minima = minima
        self._maxima = maxima
        self._extr = extr
        
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
            pos = np.greater(diff_y_min, thres)
        else:
            pos = np.greater(diff_y_max, thres)
        x = extr.x[1:-1][pos]
        y = extr.y[1:-1][pos]
        dy = extr.dy[1:-1][pos]
        xmin = []
        xmax = []        
        if len(x) > 0: 

            xmin = np.asarray(extr.x[:-2][pos]) * x.unit
            xmax = np.asarray(extr.x[2:][pos]) * x.unit

        self._bound_x = np.setxor1d(xmin, xmax)
        self._bound_y = np.interp(self._bound_x, spec.x, spec.y)
        line = Line(self.spec, x=x, y=y, xmin=xmin, xmax=xmax, dy=dy)
        self.__dict__.update(line.__dict__)


    def fit(self, x=None, y=None, dy=None, guess=None, maxfev=1000):

        if (x is None):
            x = self._fit_x
        if (y is None):
            y = self._fit_y
        if (dy is None):
            dy = self._fit_dy
        if (guess is None):
            guess = self._guess
            
        (model, param) = guess
        if (hasattr(self, '_chunk_sum')):
            chunk = self._chunk_sum
        else:
            chunk = np.full(len(self._spec.t), True)
        if (len(self._spec.x[chunk]) < len(param)):
            warnings.warn("Too few data points; skipping.")
            fit = None
        else:
            fit = model.fit(y.value, param, x=x.value, weights=1/dy.value)
            #, fit_kws={'maxfev': maxfev})
        return fit    

    def fit_auto(self, x=None, line=None, i_max=10, mode=None):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (line is not None):
            if (x is not None):
                warnings.warn("x will be used; line will be disregarded.")
            x = self.x[line]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        stop = False
        aic_old = float('inf')
        redchi_old = float('inf')        
        redchi_best = float('inf')
        i = 0
        i_best = 1
        cont_corr = 1.0
        vary = False
        self._noneb = dc(self)
        self._neb = dc(self)
        self._last_add = 0.0
        while (stop == False):
            i += 1
            noneb = dc(self._noneb)
            #try:
            #    noneb._z_fit = noneb_z_fit
            #except:
            #    pass
            fit_noneb = noneb.fit_wrap(x, vary, mode)

            neb = dc(self._neb)
            #try:
            #    neb._z_fit = neb_z_fit
            #except:
            #    pass
            fit_neb = neb.fit_wrap(x, vary, mode)

            if (noneb._redchi <= neb._redchi):
                #print("noneb")
                self_temp = dc(noneb)
                fit = fit_noneb
                #noneb_z_fit = noneb._z_fit
            else:
                #print("neb")
                self_temp = dc(neb)
                fit = fit_neb
                #neb_z_fit = neb._z_fit
                
            #self = dc(self_temp)
            self.__dict__.update(self_temp.__dict__)
            print("(%i) %3.2f;" \
                  % (i, self._redchi), end=" ", flush=True)
            stop = (self._redchi < redchi_thr) \
                   or ((self._redchi<10*redchi_thr) and (self._aic>aic_old)) \
                   or (i==i_max)
            #or (self._last_add == None) \
                   
            aic_old = self._aic
            redchi_old = self._redchi            
            if (self._redchi < redchi_best): #or 1==1):
                self_best = dc(self)
                fit_best = dc(fit)
                i_best = i

            if (stop == False):
                cont_corr = self.corr_resid(self._group, cont_corr)

        
        self = dc(self_best)
        self.__dict__.update(self_best.__dict__)
        fit = fit_best
        print("best chi-squared (%i) %3.2f, %3.2f;" \
              % (i_best, self._redchi, self._aic), end=" ", flush=True)

    def fit_list(self, list_range=None, iter_range=range(10,11), mode=None,
                 plot=True):
        if (list_range is None):
            list_range = range(len(self.t))

        self_temp = dc(self)
        x_arr = self_temp.x
        #i = 0
        group_check = 0
        group_all = None
        chunk_all = None
        for l in list_range:
            start = time.time()
            print("Redshift %i (%i/%i) (%3.4f)..." \
                  % (l+1, l+1-list_range[0], len(list_range), x_arr[l].value),
                  end=" ", flush=True)
            
            # Check if the group is new
            if (np.array_equal(self_temp.group(x=x_arr[l])[1], group_check)):
                print("same group, skipping.")
            else:
                group_check = self_temp.group(x=x_arr[l])[1]
                for i in iter_range:
                    self.fit_auto(x=x_arr[l], i_max=i, mode=mode)
                    print("time: %3.2f;" % (time.time()-start), end=" ",
                          flush=True)
                    
                    if (plot == True):
                        print("close graphs to continue.")
                        self.plot(self._group, self._chunk, mode='split')
                    else:
                        print("")
                    
    def fit_prep(self, prof='voigt', vary=False, mode=None, **kwargs):
        if (hasattr(self, '_chunk_sum')):
            where = self._chunk_sum
        else:
            where = np.full(len(self._spec.t), True)
        
        self._fit_x = self._spec.x[where]
        self._fit_y = self._spec.y[where] / self._cont.y[where]
        self._fit_dy = self._spec.dy[where] / self._cont.y[where]

        self._norm_guess = self.norm(self._group, self._chunk, vary=vary)
        if (prof == 'voigt'):
            if (mode == 'use_old'):
                if (hasattr(self, '_z_fit')):
                    diff = np.setdiff1d(self._z[self._group[1]], self._z_fit)
                    z_temp = np.append(self._z_fit.value, diff.value)
                    #N_temp = np.append(
                    #    self._N_fit, model.N_guess(self._norm, ion=self._flat.ion))
                    N_temp = np.append(self._N_fit.value, voigt_def['N'])
                    b_temp = np.append(self._b_fit.value, voigt_def['b'])
                    btur_temp = np.append(self._btur_fit.value,
                                          voigt_def['btur'])
                    z = np.sort(z_temp) * self._z_fit.unit
                    N = N_temp[np.argsort(z_temp)] * self._N_fit.unit
                    b = b_temp[np.argsort(z_temp)] * self._b_fit.unit
                    btur = btur_temp[np.argsort(z_temp)] * self._b_fit.unit
                else:
                    z = []
                    N = []
                    b = []
                    btur = []
            #"""
            else:
                try:
                    z = kwargs['z']
                except:
                    z = []
                try:
                    N = kwargs['N']
                except:
                    N = []
                try:
                    b = kwargs['b']
                except:
                    b = []
                try:
                    btur = kwargs['btur']
                except:
                    btur = []
            self._prof_guess = self.voigt(self._group, self._chunk, z=z, N=N,
                                          b=b, btur=btur)
        else:
            raise Exception("Only Voigt profile is supported.")
        self._psf = self.psf(self._group, self._chunk, self._resol)

        if (hasattr(self, '_fit') == False):
            self._fit = dc(self._spec)        
        if (hasattr(self, '_cont') == False):
            self._cont = dc(self._spec)            
        if (hasattr(self, '_resid_fit') == False):
            self._resid_fit = dc(self._spec)        
        if (hasattr(self, '_resid_cont') == False):
            self._resid_cont = dc(self._spec)            
        if (hasattr(self, '_rem') == False):
            self._rem = dc(self._spec)            

    def fit_prod(self, fit, prof='voigt'):
        if (hasattr(self, '_chunk_sum')):
            where = self._chunk_sum
        else:
            where = np.full(len(self._spec.t), True)
        yunit = self._spec.y.unit
        comp = fit.eval_components(x=self._spec.x[where].value)
        self._fit.y[where] = fit.best_fit * self._cont.y[where] 
        self._cont.y[where] = comp['cont_'] * self._cont.y[where]
        self._resid_fit.y[where] = self._spec.y[where] - self._fit.y[where]
        self._resid_cont.y[where] = self._spec.y[where] - self._cont.y[where] 
        self._rem.y[where] = self._cont.y[where] + self._resid_fit.y[where]
    #* self._spec.y[where] / self._fit.y[where] #* yunit
        
        if (prof == 'voigt'):
            z_tags = [z for z in fit.best_values if z.endswith('_z')]
            N_tags = [N for N in fit.best_values if N.endswith('_N')]
            b_tags = [b for b in fit.best_values if b.endswith('_b')]
            btur_tags = [bt for bt in fit.best_values if bt.endswith('_btur')]

            #print([fit.best_values[z] for z in z_tags])
            z_best, un = np.unique([fit.best_values[z] for z in z_tags],
                                   return_index=True)
            N_best = np.array([fit.best_values[N] for N in N_tags])[un]
            b_best = np.array([fit.best_values[b] for b in b_tags])[un]
            btur_best = np.array([fit.best_values[bt] for bt in btur_tags])[un]

            """
            self._z_fit = [fit.best_values[z] for z in np.sort(z_tags)]
            self._N_fit = [fit.best_values[N] for N in np.sort(N_tags)] \
                          * self._N_arr[0].unit
            self._b_fit = [fit.best_values[b] for b in np.sort(b_tags)] \
                          * self._b_arr[0].unit
            self._btur_fit = [fit.best_values[bt] for bt in np.sort(btur_tags)]\
                             * self._btur_arr[0].unit
            """
            self._z_fit = np.sort(z_best) * u.nm/u.nm
            self._N_fit = N_best[np.argsort(z_best)] / u.cm**2
            self._b_fit = b_best[np.argsort(z_best)] * u.km/u.s           
            self._btur_fit = btur_best[np.argsort(z_best)] * u.km/u.s
            #print(self._z_fit)
            #print(len(self._z), len(self._group[1]), len(self._z_fit))

            # When new redshift is a duplicate
            if ((hasattr(self, '_last_add')) \
                and (len(self._z_fit) < np.sum(self._group[1]))): 
                #print(self._last_add)
                self._z = np.delete(self._z, self._last_add)
                self._t.remove_row(self._last_add)
                line = self._group[0]
                self.group(line=line)
                #self._group[1] = np.delete(self._group[1], self._last_add)
            #print(len(self._z), len(self._group[1]), len(self._z_fit))
            self._z[self._group[1]] = self._z_fit
            
        else:
            raise Exception("Only Voigt profile is supported.")

        self._redchi = fit.redchi
        self._aic = fit.aic    
        
    def fit_wrap(self, x, vary=False, mode=None):
        group = self.group(x)
        chunk = self.chunk(x)
        self.fit_prep(mode=mode)
        guess = self.model()
        fit = self.fit()
        self.fit_prod(fit)
        return fit

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

    def norm(self, group, chunk, value=1.0, vary=False):
        """ Normalize continuum """

        model = Model(self._spec, line=self, group=group, chunk=chunk) 
        norm = model.norm(value, vary)
        if (hasattr(self, '_norm') == False):
            self._norm = dc(self._spec)
        self._norm.y[chunk[1]] = norm[0].eval(
            norm[1], x=self._norm.x[chunk[1]].value) \
            * self._cont.y[chunk[1]] * self._norm.y[chunk[1]].unit 

        return norm 
            
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
        
    
    def psf(self, group, chunk, resol):
        """ Model the instrumental PSF """
        
        model = Model(self._spec, line=self, group=group, chunk=chunk)
        psf = model.psf(resol)

        return psf   

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
        
    def unabs(self, group, chunk):
        """ Remove lines """

        model = Model(self._spec, line=self, group=group, chunk=chunk) 
        unabs = model.unabs()
        if (hasattr(self, '_unabs') == False):
            self._unabs = dc(self._spec)
        self._unabs.y[chunk[1]] = unabs[0].eval(
            unabs[1], x=self._unabs.x[chunk[1]].value) \
            * self._unabs.y[chunk[1]].unit

        return unabs 
    
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
