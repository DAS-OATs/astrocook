from . import Spec1D, Model
from .utils import dict_wave, voigt_def
from astropy.table import Column, Table
from astropy import units as u
from copy import deepcopy as dc
import inspect
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax
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

    def auto(self, group, chunk, unabs_guess, voigt_guess):

    def chunk(self, x=None, line=None):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x-x) == abs(self.x-x).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        iter = range(len(self._t))
        sel = self._spec.t['X'] < 0.0
        for row in self.t[self.group(line=line)[1]]:
            sel = np.logical_or(sel, np.logical_and(
                self._spec.t['X'] > row['XMIN'],
                self._spec.t['X'] < row['XMAX']))
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

    def fit(self, group, chunk, unabs_guess, voigt_guess):

        model = unabs_guess[0] * voigt_guess[0]
        param = unabs_guess[1]
        param.update(voigt_guess[1])
        if (hasattr(self, '_fit') == False):
            self._fit = dc(self._spec)
        fit = model.fit(self._spec.y[chunk[1]].value, param,
                        x=self._spec.x[chunk[1]].value,
                        weights=1/self._spec.dy[chunk[1]].value)
        self._fit.y[chunk[1]] = fit.best_fit * self._fit.y[chunk[1]].unit
        self._redchi = fit.redchi
        return fit
        
        
    def group(self, x=None, line=None):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x-x) == abs(self.x-x).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("line is too large.")
        iter = range(len(self._t))
        self._t.sort('X')
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
