from . import Line, Model
from .utils import voigt_def
from astropy import units as u
from astropy.table import Column, Table
from copy import deepcopy as dc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
import warnings
import sys

class Syst(Line):

    def __init__(self,
                 line=None,
                 spec=None,
                 x=[],
                 y=[],
                 xmin=[],
                 xmax=[],
                 dy=[],
                 ion=[],
                 yunit=None,
                 meta=None,
                 dtype=float):
        """ Constructor for the Syst class """ 
        
        # Exceptions
        #if ((ion is not []) and (line is not None) \
        #    and (len(ion) != len(line.t))):
        #    raise Exception("Ion list and line list must have the same length.")
        if ((x is []) != (y is [])):
            raise Exception("X and Y must be provided together.")
        if ((xmin is []) != (xmax is [])):
            raise Exception("XMIN and XMAX must be provided together.")
        sumlen = len(x) + len(y) + len(xmin) + len(xmax) + len(dy) + len(ion)
        if ((x != []) and (sumlen % len(x) != 0)):
            raise Exception("Data arrays must have the same length.")
        if ((xmin != []) and (sumlen % len(xmin) != 0)):
            raise Exception("Data arrays must have the same length.")

        # Warnings
        if ((spec is None) and (line is None) and (x is [])):
            warnings.warn("No spectrum, line list, or data provided.")

        # Line list
        self._line = None
        if (line is not None):
            self._line = dc(line)

        # Spectrum
        self._spec = None
        if (spec is not None):
            self._spec = dc(spec)
            
        # Ion list
        if (line is not None):
            len_ion = len(line._t)
        if (x != []):
            len_ion = len(x)

        if ((ion is []) or (ion is 'Ly_a')):
            ion = np.full(len_ion, 'Ly_a')           
        if (ion is 'CIV'):
            ion = np.stack((np.full(len_ion, 'CIV_1548'),                    
                            np.full(len_ion, 'CIV_1550'))).T
        self._ion = ion
            
        # System list
        data = ()
        if (x != []):
            col_x = Column(np.asarray(dc(x), dtype=dtype), name='X')
            col_y = Column(np.asarray(dc(y)), name='Y')
            data = (col_x, col_y)
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
            col_dy = Column(np.asarray(dc(dy)), name='DY')
            data += (col_dy,)
        if ((x != []) and (ion != [])):
            col_ion = Column(np.asarray(dc(ion)), name='ION')
            data += (col_ion,)
        if data is ():
            data = None
            
        if (meta is None):
            meta = {}

        self._t = Table(data=data, masked=True, meta=meta)
        if (y != []):
            self._t['Y'].unit = yunit
        if (dy != []):
            self._t['DY'].unit = yunit

        self._use_good = False

                
# Properties

    @property
    def ion(self):
        if self._use_good:
            ret = np.asarray(self._t['ION'][self._igood])
        else:
            ret = np.asarray(self._t['ION'])
        return ret 

    @ion.setter
    def ion(self, value):
        if self._use_good:
            self._t['ION'][self._igood] = np.asarray(value)
        else:
            self._t['ION'] = np.asarray(value)
        
    @property
    def line(self):
        return self._line

    @line.setter
    def line(self, value):
        if isinstance(value, Line):
            self._line = value
        else:
            raise Exception("Line list has a wrong format.")
    
    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, value):
        if isinstance(value, Syst):
            self._z = value
        else:
            raise Exception("Redshift list has a wrong format.")


# Methods

    def chunk(self, x=None, line=None):  # Chunk must be shifted to the system z
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x-x) == abs(self.x-x).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        try:  # When ION has different sizes in different rows
            ion = np.unique(np.asarray(np.sum(self.ion)))
        except:
            ion = np.unique(self.ion)
        n = len(ion)
        iter = range(len(self._t))
        ret = (line,)
        for p in range(n):
            sel = self._spec.t['X'] < 0.0
            spec = dc(self._spec)
            spec.to_z([ion[p]])
            for row in self.t[self.group(line=line)[1]]:
                sel = np.logical_or(sel, np.logical_and(
                    spec.t['X'] > row['XMIN'],
                    spec.t['X'] < row['XMAX']))
            ret += (sel,)
        return ret

    def create_z(self):

        # Redshift list
        self._z = dc(self._line)
        self._z.to_z(self._ion)
        self._z.t.add_column(Column(np.asarray(dc(self._ion)), name='ION'))

        if self._use_good:
            self._z.ion = np.asarray(self._z.t['ION'][self._igood])
        else:
            self._z.ion = np.asarray(self._z.t['ION'])
            
    def match_z(self, ztol=1e-4):
        """ Match redshifts in a list, to define systems """

        if (hasattr(self, '_z') == False):
            raise Exception("Redshift table must be created before matching.")
        
        
        # Flatten arrays
        # N.B. Y and DY don't need flattening, as they have only one entry
        # per row. They will be repeated to have the shape of the other arrays.
        z = np.ravel(self._z.x)
        y = np.repeat(self._z.y, self._z.x.shape[1])        
        zmin = np.ravel(self._z.xmin)
        zmax = np.ravel(self._z.xmax)
        dy = np.repeat(self._z.dy, self._z.x.shape[1])
        ion = np.ravel(self._z.ion)

        # Sort arrays
        argsort = np.argsort(z, kind='mergesort')
        z_sort = np.asarray(z[argsort])
        y_sort = np.asarray(y[argsort])      
        zmin_sort = np.asarray(zmin[argsort])
        zmax_sort = np.asarray(zmax[argsort])
        dy_sort = np.asarray(dy[argsort])
        ion_sort = np.asarray(ion[argsort])
        
        # Find coincidences
        coinc = np.isclose(z_sort[1:], z_sort[:-1], atol=ztol) 
        coupl = np.core.defchararray.not_equal(ion_sort[1:], ion_sort[:-1])
        coinc = np.logical_and(coinc, coupl)

        new = True
        z_coinc = []
        y_coinc = []
        zmin_coinc = []        
        zmax_coinc = []
        ion_coinc = []
        dy_coinc = []
        for i in range(len(z_sort) - 2):
            if ((new == True) and (coinc[i] == True)):
                new = False
                z_coinc_row = z_sort[i]
                y_coinc_row = (y_sort[i],)
                zmin_coinc_row = zmin_sort[i]
                zmax_coinc_row = zmax_sort[i]                
                dy_coinc_row = (dy_sort[i],)
                ion_coinc_row = (ion_sort[i],)
            if (coinc[i] == True):
                z_coinc_row = np.append(z_coinc_row, z_sort[i+1])
                y_coinc_row = y_coinc_row + (y_sort[i+1],)
                zmin_coinc_row = np.append(zmin_coinc_row, zmin_sort[i+1])
                zmax_coinc_row = np.append(zmax_coinc_row, zmax_sort[i+1])
                dy_coinc_row = dy_coinc_row + (dy_sort[i+1],)
                ion_coinc_row = ion_coinc_row + (ion_sort[i+1],)
            if (coinc[i+1] == False):
                if (new == False):
                    z_coinc.append(np.mean(z_coinc_row))
                    y_coinc.append(y_coinc_row)
                    zmin_coinc.append(np.mean(zmin_coinc_row))
                    zmax_coinc.append(np.mean(zmax_coinc_row))
                    dy_coinc.append(dy_coinc_row)
                    ion_coinc.append(ion_coinc_row)
                new = True
        syst = Syst(self.line, self.spec, x=z_coinc, y=y_coinc,
                    xmin=zmin_coinc, xmax=zmax_coinc, dy=dy_coinc,
                    ion=ion_coinc, yunit=y.unit)
        self.__dict__.update(syst.__dict__)

    def fit(self, group, chunk, unabs_guess, voigt_guess):


        for c in range(1, len(chunk)):
            if (c == 1):
                model = unabs_guess[2*(c-1)] * voigt_guess[2*(c-1)]
                param = unabs_guess[2*c-1]
                chunk_sum = chunk[c]
            else:
                model += unabs_guess[2*(c-1)] * voigt_guess[2*(c-1)]
                param.update(unabs_guess[2*c-1])
                chunk_sum += chunk[c]
            param.update(voigt_guess[2*c-1])
        expr_dict = voigt_guess[-1]
        for k in expr_dict:
            param[k].set(expr=expr_dict[k])

        if (hasattr(self, '_fit') == False):
            self._fit = dc(self._spec)

        fit = model.fit(self._spec.y[chunk_sum].value, param,
                        x=self._spec.x[chunk_sum].value,
                        weights=1/self._spec.dy[chunk_sum].value)
        
        self._fit.y[chunk_sum] = fit.best_fit * self._fit.y[chunk_sum].unit
        self._redchi = fit.redchi
        
        return fit    

    def flatten_z(self):
        """ Create a flattened version of the system, with different entries
        for each ion """

        yunit = self._z.dy.unit
        first = True
        for r in self.t:
            for i in range(len(r['Y'])):
                if (first == True):
                    self._flat = Syst(x=[r['X']], y=[r['Y'][i]],
                                      xmin=[r['XMIN']], xmax=[r['XMAX']],
                                      dy=[r['DY'][i]], ion=[r['ION'][i]],
                                      yunit=yunit)
                    first = False
                else:
                    self._flat.t.add_row([r['X'], r['Y'][i], r['XMIN'],
                                          r['XMAX'], r['DY'][i], r['ION'][i]])
                    
    def plot(self, group=None, chunk=None, figsize=(6,6), split=False,
             block=True, **kwargs):
        ion = np.unique(self._flat.ion)
        n = len(ion)
        z = self.x
        if (chunk is not None):
            zmin = np.min(self.xmin[group[1]])
            zmax = np.max(self.xmax[group[1]])
        else:
            zmin = np.min(self.xmin)
            zmax = np.max(self.xmax)
        if split == True:
            row = min(n,4)
            col = int(np.ceil(n/4))
            figsize = (col*6, n*3.5)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            grid = gs(row,col)            
            for p in range(n):
                ax = fig.add_subplot(grid[p%4, int(np.floor(p/4))])
                ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
                spec = dc(self._spec)
                line = dc(self._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
                ax.plot(spec.x, spec.y, c='b')
                if (hasattr(self, '_unabs')):
                    unabs = dc(self._unabs)
                    unabs.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        ax.plot(unabs.x[chunk[c]], unabs.y[chunk[c]], c='y',
                                linestyle=':')
                if (hasattr(self, '_voigt')):
                    voigt = dc(self._voigt)
                    voigt.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        ax.plot(voigt.x[chunk[c]], voigt.y[chunk[c]], c='g',
                                linestyle=':')
                if (hasattr(self, '_fit')):
                    fit = dc(self._fit)
                    fit.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        ax.plot(fit.x[chunk[c]], fit.y[chunk[c]], c='g')
                ax.scatter(line.x, line.y, c='b')
                for comp in z:
                    ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='black')
                if ((p+1) % row != 0):
                    ax.set_xticks([], [])
                else:
                    ax.set_xlabel("Redshift")
                ax.set_xlim(zmin, zmax)
                ax.text(0.5, 0.92, ion[p], horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes,
                        fontsize=12)
            if (hasattr(self, '_fit')):
                fig.suptitle("Reduced chi-squared: %3.2f" % (self._redchi),
                             fontsize=10)
        else:
            grid = gs(1,1)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            ax = fig.add_subplot(grid[:,:])
            ax.set_xlabel("Redshift")
            ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
            for p in range(n):
                spec = dc(self._spec)
                line = dc(self._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
                ax.set_xlim(zmin, zmax)
                ax.plot(spec.x, spec.y)
                ax.scatter(line.x, line.y)
            text = ', '.join(str(p) for p in ion)
            for comp in z:
                ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='black')

            ax.text(0.5, 0.92, text, horizontalalignment="center",
                    verticalalignment="center", transform=ax.transAxes,
                    fontsize=12)
           
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
        grid.update(hspace=0.0)
        if block is False:
            plt.ion()
            plt.draw()
        else:
            plt.show()

    def unabs(self, group, chunk):
        """ Remove lines """

        model = Model(self._spec, syst=self, group=group, chunk=chunk) 
        unabs = model.unabs()
        if (hasattr(self, '_unabs') == False):
            self._unabs = dc(self._spec)
        for c in range(1, len(chunk)):
            self._unabs.y[chunk[c]] = unabs[2*(c-1)].eval(
                unabs[2*c-1], x=self._unabs.x[chunk[c]].value) \
                * self._unabs.y[chunk[c]].unit
    
        return unabs

    def voigt(self, group, chunk, z=[], N=[], b=[], btur=[]):

        sumlen = len(z) + len(N) + len(b) + len(btur)
        if ((z != []) and (sumlen % len(z) != 0)):
            raise Exception("Parameter arrays must have the same length.")

        model = Model(self._spec, syst=self, group=group, chunk=chunk)
        if (z == []):
            z = self.x[group[1]]
            N = model.N_guess(self._unabs, ion=self._flat.ion)
            b = np.full(len(self.x[group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(self.x[group[1]]), voigt_def['btur']) \
                   * u.km / u.s

        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)

        ion = np.unique(self._flat.ion)
        voigt = model.voigt(z, N, b, btur, ion)
        for c in range(1, len(chunk)):
            self._voigt.y[chunk[c]] = voigt[2*(c-1)].eval(
                voigt[2*c-1], x=self._voigt.x[chunk[c]].value) \
                * self._unabs.y[chunk[c]]
        return voigt
   
