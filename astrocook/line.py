from . import Spec1D, Model
from .utils import convolve, dict_wave, redchi_thr, voigt_def
from astropy.table import Column, Table
from astropy import units as u
from copy import deepcopy as dc
import inspect
from lmfit import CompositeModel as lmc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve
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
        if ((x is []) != (y is [])):
            raise Exception("X and Y must be provided together.")
        if ((xmin is []) != (xmax is [])):
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

    def add_min_resid(self, group):
        """ Add a new line at the minimum residual """

        xmin = np.min(self.xmin[group[1]])
        xmax = np.max(self.xmax[group[1]])
        where = np.where(np.logical_and(self._resid.x >= xmin,
                                        self._resid.x <= xmax))
        resid_norm = self._resid.y.value * 0.0
        resid_norm[where] = self._resid.y[where]/self._resid.dy[where]
        x = self._resid.x[np.argmin(resid_norm)]
        y = np.interp(x.value, self._spec.x, self._spec.y) 
        dy = np.interp(x.value, self._spec.x, self._spec.dy)
        self._t.add_row([x, y, xmin, xmax, dy])

    def auto(self, x=None, line=None):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x.value-x.value) \
                            == abs(self.x.value-x.value).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        stop = False
        aic_old = float('inf')
        redchi_old = float('inf')        
        redchi_min = float('inf')
        i = 0
        print("")
        print(" Chi-squared:", end=" ", flush=True)
        while (stop == False):
            i += 1
            group = self.group(line=line)            
            chunk = self.chunk(line=line)
            unabs_guess = self.unabs(group, chunk)
            voigt_guess = self.voigt(group, chunk)
            psf = self.psf(group, chunk, self._resol)
            fit = self.fit(group, chunk, unabs_guess, voigt_guess, psf)
            self._resid = dc(self._spec)
            self._resid.y = self._spec.y - self._fit.y
            
            print("(%i) %3.2f, %3.2f;" % (i, self._redchi, self._aic), end=" ",
                  flush=True)
            stop = (self._redchi < redchi_thr) \
                   or ((self._redchi<10*redchi_thr) and (self._aic>aic_old)) \
                   or (self._redchi>1e4) \
                   or (i==10)

            aic_old = self._aic
            redchi_old = self._redchi            
            if (self._redchi < redchi_min):
                redchi_min = self._redchi
                group_best = group
                chunk_best = chunk
                unabs_best = unabs_guess
                voigt_best = voigt_guess
                psf = psf
                i_best = i
                t_best = dc(self._t)
            if (stop == False):
                self.add_min_resid(group)
                
            
        self._t = dc(t_best)       
        fit = self.fit(group_best, chunk_best, unabs_best, voigt_best, psf)
        print("best: (%i) %3.2f, %3.2f;" % (i_best, self._redchi, self._aic),
              end=" ", flush=True)
        
        return group_best, chunk_best

    def chunk(self, x=None, line=None):
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
        for row in self.t[self.group(line=line)[1]]:
            sel = np.logical_or(sel, np.logical_and(
                self._spec.t['X'] >= row['XMIN'],
                self._spec.t['X'] <= row['XMAX']))
        if (np.sum(sel) % 2 == 0):
            sel[np.argmax(sel)] = 0
        return (line, sel)

    def find(self, mode='abs', diff='max', kappa=3.0, hwidth=2):
        """ Find lines in a spectrum """

        if (self._spec is None):
            raise Exception("Spectrum not provided.")
        
        spec = self._spec
        
        # Find extrema
        minima, maxima, extr = spec.find_extrema()

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

            # Compute the boundaries for line fitting    
            window = extr.rolling_window(hwidth * 2 + 1)
            bound = np.full(hwidth + 1, np.amin(spec.x))
            if mode is 'em':
                app = window.x[np.arange(len(window.x)),
                               np.argmin(window.y, 1)]
            else:
                app = window.x[np.arange(len(window.x)),
                               np.argmax(window.y, 1)]
            bound = np.append(bound, app)
            bound = np.append(bound, np.full(hwidth + 1, np.amax(spec.x)))
            
            xmin = np.asarray(bound[:-hwidth * 2][pos]) * x.unit 
            xmax = np.asarray(bound[hwidth * 2:][pos]) * x.unit

        line = Line(self.spec, x=x, y=y, xmin=xmin, xmax=xmax, dy=dy)
        self.__dict__.update(line.__dict__)

    def fit(self, group, chunk, unabs_guess, voigt_guess, psf, maxfev=300):

        model = unabs_guess[0] * voigt_guess[0]
        param = unabs_guess[1]
        conv_model = lmc(model, psf[0], convolve)
        param.update(voigt_guess[1])
        param.update(psf[1])
        if (hasattr(self, '_fit') == False):
            self._fit = dc(self._spec)
        if (hasattr(self, '_cont') == False):
            self._cont = dc(self._spec)
        if (len(self._spec.x[chunk[1]]) < len(param)):
            warnings.warn("Too few data points; skipping.")
            fit = None
        else:
            fit = conv_model.fit(self._spec.y[chunk[1]].value, param,
                                 x=self._spec.x[chunk[1]].value,
                                 fit_kws={'maxfev': maxfev},
                                 weights=1/self._spec.dy[chunk[1]].value)
            cont = fit.eval_components(x=self._spec.x[chunk[1]].value)
            self._cont.y[chunk[1]] = cont['cont1_'] \
                                     * self._cont.y[chunk[1]].unit
            self._fit.y[chunk[1]] = fit.best_fit * self._fit.y[chunk[1]].unit
            self._redchi = fit.redchi
            self._aic = fit.aic

        return fit
        
    def group(self, x=None, line=None):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x.value-x.value) \
                            == abs(self.x.value-x.value).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("line is too large.")
        iter = range(len(self._t))
        self._t.sort('X')  # This gives a warning on a masked table
        #self._t.group_by('X')
        groups = np.append([0], np.cumsum(self._t['XMAX'][:-1] <
                                          self._t['XMIN'][1:])) 
        sel = np.array([groups[l] == groups[line] for l in iter])
        return (line, sel)
        
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

    def plot(self, group=None, chunk=None, figsize=(10,4), block=True,
             **kwargs):
        spec = self._spec
        fig = plt.figure(figsize=figsize)
        fig.canvas.set_window_title("Lines")
        grid = gs(1, 1)
        ax = fig.add_subplot(grid[:, :])
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
        ax.plot(spec.x, spec.y, c='b')
        if (chunk is not None):
            ax.plot(spec.x[chunk[1]], spec.y[chunk[1]], c='r')
            if (hasattr(self, '_unabs')):
                where = np.where(self._spec.y != self._unabs.y)
                ax.plot(self._unabs.x[where], self._unabs.y[where], c='y',
                        linestyle=':')
            if (hasattr(self, '_voigt')):
                where = np.where(self._spec.y != self._voigt.y)
                ax.plot(self._voigt.x[where], self._voigt.y[where], c='g',
                        linestyle=':')
            if (hasattr(self, '_cont')):
                where = np.where(self._spec.y != self._cont.y)
                ax.plot(self._cont.x[where], self._cont.y[where], c='y')
            if (hasattr(self, '_fit')):
                where = np.where(self._spec.y != self._fit.y)
                ax.plot(self._fit.x[where], self._fit.y[where], c='g')
        ax.scatter(self.x, self.y, c='b')
        if (group is not None):
            ax.scatter(self.x[group[0]], self.y[group[0]], c='r')        
            ax.scatter(self.x[group[1]], self.y[group[1]], c='g')        
        if (hasattr(self, '_fit')):
            fig.suptitle("Reduced chi-squared: %3.2f" % (self._redchi),
                         fontsize=10)
        ax.set_xlabel("Wavelength [" + str(self._spec.x.unit) + "]")
        ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")

        
        if block is False:
            plt.ion()
            plt.draw()
        else:
            plt.show()

    def psf(self, group, chunk, resol):
        """ Model the instrumental PSF """
        
        model = Model(self._spec, line=self, group=group, chunk=chunk)
        psf = model.psf(resol)

        return psf   

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
            N = model.N_guess(self._unabs)
            b = np.full(len(self.x[group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(self.x[group[1]]), voigt_def['btur']) \
                   * u.km / u.s
        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)
        voigt = model.voigt(z, N, b, btur, 'Ly_a')
        self._voigt.y[chunk[1]] = voigt[0].eval(
            voigt[1], x=self._voigt.x[chunk[1]].value) \
            * self._unabs.y[chunk[1]]

        return voigt
