from . import Line, Model
from .utils import convolve, convolve2, dict_doubl, dict_wave, redchi_thr, \
    voigt_def
from .model import voigt_params
from astropy import units as u
from astropy.io import fits as fits
from astropy.table import Column, Table
from copy import deepcopy as dc
from lmfit import CompositeModel as lmc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
import warnings
import sys
import time

class System(Line):

    def __init__(self,
                 line=None,
                 spec=None,
                 x=[],
                 y=[],
                 xmin=[],
                 xmax=[],
                 dy=[],
                 ion=[],
                 doubl=[],
                 N=[],
                 b=[],
                 btur=[],
                 yunit=None,
                 meta=None,
                 dtype=float):
        """ Constructor for the System class """ 
        
        # Exceptions
        if ((x is []) != (y is [])):
            raise Exception("X and Y must be provided together.")
        if ((xmin is []) != (xmax is [])):
            raise Exception("XMIN and XMAX must be provided together.")
        sumlen = len(x) + len(y) + len(xmin) + len(xmax) + len(dy) + len(ion) \
                 + len(doubl) 
        if ((x != []) and (sumlen % len(x) != 0)):
            raise Exception("Data arrays must have the same length.")
        if ((xmin != []) and (sumlen % len(xmin) != 0)):
            raise Exception("Data arrays must have the same length.")
        sumlen_voigt = len(x) + len(N) + len(b) + len(btur)
        if ((N != []) and (sumlen_voigt % len(x) != 0)):
            raise Exception("Voigt parameters must have the same length of X.")

        
        # Warnings
        if ((spec is None) and (line is None) and (x is [])):
            warnings.warn("No spectrum, line list, or data provided.")

        # Line list
        self._line = None
        if (line is not None):
            self._line = dc(line)
            if (hasattr(line, '_cont')):
                self._precont = dc(line._precont)
                self._cont = dc(line._cont)
            if (hasattr(line, '_minima')):
                self._minima = dc(line._minima)
            if (hasattr(line, '_maxima')):
                self._maxima = dc(line._maxima)            

        # Spectrum
        self._spec = None
        if (spec is not None):
            """
            if (hasattr(spec, '_orig')):
                self._spec = dc(spec._orig)
            else:
                self._spec = dc(spec)
            """
            self._spec = dc(spec)
            
        # Ion list
        if (line is not None):
            len_ion = len(line._t)
        if (x != []):
            len_ion = len(x)

        if ((ion == []) and (doubl != [])):
            ion = [dict_doubl[doubl]]
            for i in range(1, len_ion):
                ion = np.append(ion, [dict_doubl[doubl]], 0)
            #ion = np.stack((np.full(len_ion, dict_doubl[doubl][0]),
            #                np.full(len_ion, dict_doubl[doubl][1]))).T
                           
        if ((ion is []) or (ion is 'Ly_a')):
            ion = np.full(len_ion, 'Ly_a')           

        if (N == []):
            N = np.full(len_ion, float('nan')) 
            b = np.full(len_ion, float('nan')) 
            btur = np.full(len_ion, float('nan'))
        dN = np.full(len_ion, float('nan')) 
        db = np.full(len_ion, float('nan')) 
        dbtur = np.full(len_ion, float('nan'))
        
        self._ion = ion

        # System list
        data = ()
        if (x != []):
            col_x = Column(np.asarray(dc(x), dtype=dtype), name='X')
            col_y = Column(np.asarray(dc(y)), name='Y', unit=yunit)
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
            col_dy = Column(np.asarray(dc(dy)), name='DY', unit=yunit)
            data += (col_dy,)
        if ((x != []) and (ion != [])):
            col_ion = Column(np.asarray(dc(ion)), name='ION')
            data += (col_ion,)
        if ((x != []) and (ion != [])):
            col_N = Column(np.asarray(dc(N)), name='N', unit=1/u.cm**2)
            col_b = Column(np.asarray(dc(b)), name='B', unit=u.km/u.s)
            col_btur = Column(np.asarray(dc(btur)), name='BTUR',
                              unit=u.km/u.s)
            col_dN = Column(np.asarray(dc(dN)), name='DN', unit=1/u.cm**2)
            col_db = Column(np.asarray(dc(db)), name='DB', unit=u.km/u.s)
            col_dbtur = Column(np.asarray(dc(dbtur)), name='DBTUR',
                               unit=u.km/u.s)
            data += (col_N, col_b, col_btur, col_dN, col_db, col_dbtur)

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
    def linez(self):
        return self._linez

    @linez.setter
    def linez(self, value):
        if isinstance(value, System):
            self._linez = value
        else:
            raise Exception("Redshift list has a wrong format.")

    @property
    def t(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t        

# Methods

    def add_comp(self, cont_corr):
        """ Add a component to a line group

        The component is added at the position of the strongest negative 
        residual.
        """

        where = self._chunk_sum
        resid_norm = np.full(len(self._resid_fit.y.value), 
                             np.max(self._resid_fit.y[where]\
                                    /self._resid_fit.dy[where]) * 10)
        resid_norm[where] = self._resid_fit.y[where]/self._resid_fit.dy[where]
        x = self._resid_fit.x[where][np.argmin(resid_norm[where])]
        y = np.interp(x.value, self._spec.x, self._spec.y) * self._spec.yunit
        dy = np.interp(x.value, self._spec.x, self._spec.dy) * self._spec.yunit

        ion_arr = np.unique(self._flat.ion)
        n = len(ion_arr)
        z_arr = np.empty(n)        
        z_cen = np.empty(n)        

        for p in range(n):
            z_arr[p] = x / dict_wave[ion_arr[p]] - 1.0
            tab = Table(self.t[self._group[1]][ion_arr[p] in 'ION'])
            z_cen[p] = tab['X']

        where = (abs(z_arr-z_cen) == abs(z_arr-z_cen).min())
        z_ion = z_arr[where][0]

        ion_where = (abs(z_ion-self.x) == abs(z_ion-self.x).min())
        ion = self.ion[ion_where][0]
        zmin_ion = self.xmin[ion_where][0] 
        zmax_ion = self.xmax[ion_where][0]

        size = np.size(ion)
        for i in range(size):
            spec = dc(self._spec)
            if (size > 1):
                spec.to_z([ion[i]])
            else:
                spec.to_z([ion])
            #y_ion[i] = np.interp(z_ion, spec.x, spec.y)
            #dy_ion[i] = np.interp(z_ion, spec.x, spec.dy)
            if (i == 0):
                y_ion = (np.interp(z_ion, spec.x, spec.y),)
                dy_ion = (np.interp(z_ion, spec.x, spec.dy),)
            else:
                y_ion += (np.interp(z_ion, spec.x, spec.y),)
                dy_ion += (np.interp(z_ion, spec.x, spec.dy),)
                
        self._z_add = z_ion
        self._t.add_row([z_ion, y_ion, zmin_ion, zmax_ion, dy_ion, 
                        ion, float('nan'), float('nan'), float('nan'),
                        float('nan'), float('nan'), float('nan')])
        self._t.sort('X')  # This gives an annoying warning
        #print(self._z)
        #self._z = np.unique(np.append(self._z, z_ion)) 
        #print(self._z)
        #self._z = np.append(self._z, z_ion) 
        #print(self._z)
        #self._z.sort()
        #self._z *= z_ion.unit
        self._z = self._t['X'] * u.nm/u.nm
        self._last_add = np.where(self._z == z_ion)[0][0]
                
        
    def corr_resid(self, cont_corr):
        """ Add a new line at the minimum residual """

        neb = 'neb'
        
        where = self._chunk_sum
        resid_norm = np.full(len(self._resid_fit.y.value), 
                             np.max(self._resid_fit.y[where]\
                                    /self._resid_fit.dy[where]) * 10)
        resid_norm[where] = self._resid_fit.y[where]/self._resid_fit.dy[where]
        x = self._resid_fit.x[where][np.argmin(resid_norm[where])]
        y = np.interp(x.value, self._spec.x, self._spec.y) * self._spec.yunit
        dy = np.interp(x.value, self._spec.x, self._spec.dy) * self._spec.yunit

        ion_arr = np.unique(self._flat.ion)
        n = len(ion_arr)
        z_arr = np.empty(n)        
        z_cen = np.empty(n)        

        #xmin_ion = float('inf') * u.nm
        #xmax_ion = 0 * u.nm

        for p in range(n):
            z_arr[p] = x / dict_wave[ion_arr[p]] - 1.0
            tab = Table(self.t[self._group[1]][ion_arr[p] in 'ION'])
            z_cen[p] = tab['X']

        where = (abs(z_arr-z_cen) == abs(z_arr-z_cen).min())
        z_ion = z_arr[where][0]

        ion_where = (abs(z_ion-self.x) == abs(z_ion-self.x).min())
        ion = self.ion[ion_where][0]
        zmin_ion = self.xmin[ion_where][0] 
        zmax_ion = self.xmax[ion_where][0]
        #y_ion = np.empty(len(ion))
        #dy_ion = np.empty(len(ion))        

        size = np.size(ion)
        for i in range(size):
            spec = dc(self._spec)
            if (size > 1):
                spec.to_z([ion[i]])
            else:
                spec.to_z([ion])
            #y_ion[i] = np.interp(z_ion, spec.x, spec.y)
            #dy_ion[i] = np.interp(z_ion, spec.x, spec.dy)
            if (i == 0):
                y_ion = (np.interp(z_ion, spec.x, spec.y),)
                dy_ion = (np.interp(z_ion, spec.x, spec.dy),)
            else:
                y_ion += (np.interp(z_ion, spec.x, spec.y),)
                dy_ion += (np.interp(z_ion, spec.x, spec.dy),)
                
        
        z_neb = x / dict_wave[neb] - 1.0
        if self._line != None:
            line = dc(self._line)
            line.to_z([neb])
        else:
            line = self
        neb_where = abs(z_neb-line.x) == abs(z_neb-line.x).min()
        zmin_neb = line.xmin[neb_where][0] 
        zmax_neb = line.xmax[neb_where][0]
        spec = dc(self._spec)
        spec.to_z([neb])
        y_neb = np.interp(x.value, self._spec.x, self._spec.y) \
                * self._spec.yunit
        dy_neb = np.interp(x.value, self._spec.x, self._spec.dy) \
                 * self._spec.yunit

        self._z_add = z_ion
        
        self._noneb = dc(self)
        #if (z_ion not in self._noneb.x):
        self._noneb.t.add_row([z_ion, y_ion, zmin_ion, zmax_ion, dy_ion, 
                               ion, float('nan'), float('nan'), float('nan'),
                               float('nan'), float('nan'), float('nan')])
        self._noneb.t.sort('X')  # This gives an annoying warning
        self._noneb._z = np.unique(np.append(self._z.value, z_ion)) 
        self._noneb._z = np.append(self._z.value, z_ion) 
        self._noneb._z.sort()
        self._noneb._z *= self._z.unit
        self._noneb._last_add = np.where(self._noneb._z == z_ion)[0][0]
        #else:
        #    self._noneb._last_add = None
            
        self._neb = dc(self)
        #if (z_ion not in self._noneb.x):
        self._neb.flatten_z()
        self._neb._flat.t.add_row([z_neb, y_neb, zmin_neb, zmax_neb, dy_neb,
                                   neb, float('nan'), float('nan'),
                                   float('nan'), float('nan'), float('nan'),
                                   float('nan')])
        self._neb.deflatten_z()
        self._neb.t.sort('X')  # This gives an annoying warning
        self._neb._z = np.unique(np.append(self._z.value, z_neb.value))
        #self._neb._z = np.append(self._z.value, z_neb.value)
        self._neb._z.sort()
        self._neb._z *= self._z.unit
        self._neb._last_add = np.where(self._neb._z == z_neb)[0][0]
        #else:
        #    self._neb._last_add = None
       
        return cont_corr
        
    def chunk(self, x=None, line=None, single=False):  # Chunk must be shifted to the system z
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            #line = np.where(abs(self.x-x.value) \
            #                == abs(self.x-x.value).min())[0][0]
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
            for row in self.t[self.group(line=line, single=single)[1]]:
                sel = np.logical_or(sel, np.logical_and(
                    spec.t['X'] >= row['XMIN'],
                    spec.t['X'] <= row['XMAX']))
            if (np.sum(sel) % 2 == 0):
                sel[np.argmax(sel)] = 0
            ret += (sel,)

        self._chunk = ret
        for c in range(1, len(ret)):
            if (c == 1):
                self._chunk_sum = dc(ret[c])
            else:
                self._chunk_sum += ret[c]

        return ret

    def create_z(self):

        # Redshift list
        self._linez = dc(self._line)
        self._linez.to_z(self._ion)
        self._linez.t.add_column(Column(np.asarray(dc(self._ion)), name='ION'))

        if self._use_good:
            self._linez.ion = np.asarray(self._linez.t['ION'][self._igood])
        else:
            self._linez.ion = np.asarray(self._linez.t['ION'])

    def deflatten_z(self):
        """ Create a non-flattened version of the system from a flattened one"""

        try:
            yunit = self.y.unit
        except:
            yunit = self._line.y.unit
        self._flat.t.sort('X')  # This gives an annoying warning        
        
        first = True
        z_deflat = []
        y_deflat = []
        zmin_deflat = []        
        zmax_deflat = []
        ion_deflat = []
        dy_deflat = []
        N_deflat = []
        b_deflat = []
        btur_deflat = []        
        dN_deflat = []
        db_deflat = []
        dbtur_deflat = []        
        
        z = 0
        end_row = False
        for l in range(len(self._flat.t)):
            if (np.isclose(self._flat.x[l], z, rtol=1e-6) == False):
                if (end_row == True):
                    z_deflat.append(z)
                    y_deflat.append(y)
                    zmin_deflat.append(zmin)
                    zmax_deflat.append(zmax)
                    dy_deflat.append(dy)
                    ion_deflat.append(ion)
                    N_deflat.append(N)
                    b_deflat.append(b)
                    btur_deflat.append(btur)
                    dN_deflat.append(dN)
                    db_deflat.append(db)
                    dbtur_deflat.append(dbtur)
                    end_row = False
                y = (self._flat.y[l].value,)
                zmin = self._flat.xmin[l]
                zmax = self._flat.xmax[l]
                dy = (self._flat.dy[l].value,)                
                ion = (self._flat.ion[l],)
                N = self._flat.t['N'][l]
                b = self._flat.t['B'][l]
                btur = self._flat.t['BTUR'][l]
                dN = self._flat.t['DN'][l]
                db = self._flat.t['DB'][l]
                dbtur = self._flat.t['DBTUR'][l]
            else:
                z = self._flat.x[l]
                y = y + (self._flat.y[l].value,)
                dy = dy + (self._flat.dy[l].value,)
                ion = ion + (self._flat.ion[l],)
            z = self._flat.x[l]
            end_row = True
        z_deflat.append(z)
        y_deflat.append(y)
        zmin_deflat.append(zmin)
        zmax_deflat.append(zmax)
        dy_deflat.append(dy)
        ion_deflat.append(ion)
        N_deflat.append(N)
        b_deflat.append(b)
        btur_deflat.append(btur)
        dN_deflat.append(dN)
        db_deflat.append(db)
        dbtur_deflat.append(dbtur)

        syst = System(self.line, self.spec, x=z_deflat, y=y_deflat,
                    xmin=zmin_deflat, xmax=zmax_deflat, dy=dy_deflat,
                    ion=ion_deflat, N=N_deflat, b=b_deflat, btur=btur_deflat,
                    yunit=yunit)
        self.__dict__.update(syst.__dict__)

    def find(self, zstart=None, zend=None, ztol=1e-4, match=True):
        self.create_z()
        if (match == True):
            self.match_z(zstart, zend, ztol)
        else:
            self._z = u.Quantity(np.array(self._linez.t['X']))
            syst = System(self.line, self.spec, x=self._linez.t['X'],
                          y=self._linez.t['Y'], xmin=self._linez.t['XMIN'],
                          xmax=self._linez.t['XMAX'], dy=self._linez.t['DY'],
                          ion=self._linez.t['ION'],
                          yunit=self._linez.t['Y'].unit)
            self.__dict__.update(syst.__dict__)
        self.flatten_z()

        
    def fit(self, x=None, y=None, dy=None, guess=None, chunk=None, maxfev=1000):
        """ Fit a group of lines """
        
        if (guess is None):
            guess = self._guess
        if (chunk is None):
            chunk_sum = self._chunk_sum
        cont = self._cont.y[chunk_sum]
        slope = cont / np.mean(cont)
        x = self._spec.x[chunk_sum]
        y = self._spec.y[chunk_sum] / cont
        dy = self._spec.dy[chunk_sum] / cont
        
        # Create the model
        (model, param) = guess

        # Check if there are enough data points
        if (len(x) < len(param)):
            warnings.warn("Too few data points; skipping.")
            fit = None

        # Otherwise fit
        else:
            #print("before fit")
            #param.pretty_print()
            #print(param)
            fit = model.fit(y.value, param, x=x.value, weights=1/dy.value)#, #)
                            #fit_kws={'maxfev': maxfev})
            #print("after fit")                
            #fit.params.pretty_print()
            #print(fit.params)
        return fit

    def fit_add(self, x=None, line=None, i_max=10, mode=None, **kwargs):
        """ Fit a group of lines 
        
        Given a line, the whole group of adjacient lines is fitted, adding
        new lines if needed.
        """

        # Check input
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (line is not None):
            if (x is not None):
                warnings.warn("x will be used; line will be disregarded.")
            x = self.x[line]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")

        # Initialize 
        stop = False
        aic_old = float('inf')
        redchi_old = float('inf')        
        redchi_best = float('inf')
        i = 0
        i_best = 1
        cont_corr = 1.0
        vary = False  # To change the continuum
        self._last_add = 0.0
        while (stop == False):
            i += 1

            # Fit the group
            group = self.group(x)
            chunk = self.chunk(x)

            # Prepare the parameters
            if i == 1:
                start = voigt_params(self, **kwargs)
            else:
                start = {'z': [], 'N': [], 'b': [], 'btur': []}
            self.fit_prep(mode=mode, vary=vary, **start)

            # Create a model
            guess = self.model()

            # Fit the model
            fit = self.fit()

            if (fit == None):
                stop = True
            else:

                # Evaluate the fit
                stop = self.fit_eval(fit)

                # Save the products
                self.fit_prod(fit)

            
                # Print the current fit result
                print("(%i) %3.2f;" \
                      % (i, self._redchi), end=" ", flush=True)

                # If the current fit is the best, save it
                if (self._redchi < redchi_best): 
                    self_best = dc(self)
                    fit_best = dc(fit)
                    i_best = i
                    redchi_best = self._redchi
                    
                #"""
                # Check if iteration must stop
                cond = []
                cond.append(self._redchi < redchi_thr)
                cond.append((self._redchi<10*redchi_thr) \
                            and (self._aic>aic_old))
                cond.append(i==i_max)
                stop = np.any(cond)
                aic_old = self._aic
                redchi_old = self._redchi
                #"""

                # If not, add a line to the group
                if (stop == False):
                    cont_corr = self.add_comp(cont_corr)

        # Export best fit
        self = dc(self_best)
        self.__dict__.update(self_best.__dict__)
        fit = fit_best
        print("best chi-squared (%i) %3.2f, %3.2f; " \
              % (i_best, redchi_best, self._aic), end=" ", flush=True)

    def fit_eval(self, fit):
        """ Evaluate the improvement obtained by a fit """

        stop = False
        
        chunk = self._chunk
        chunk_sum = self._chunk_sum
        x = self._spec.x[chunk_sum]
        y = self._spec.y[chunk_sum]        
        dy = self._spec.dy[chunk_sum]

        # Division by comp['cont_'] is needed because at this point the
        # continuum has already been updated 
        comp = fit.eval_components(x=self._spec.x[chunk_sum].value)
        cont = dc(self._cont.y)
        y_fit = fit.best_fit * cont[chunk_sum]
        cont[chunk_sum] = cont[chunk_sum] * comp['cont_']
        #plt.plot(x, y)
        #plt.plot(x, y_fit)
        redchi = np.array([self.redchi(y, dy, y_fit, len(y)-fit.nvarys)])

        for c in range(1, len(chunk)):
            y_temp = dc(self._spec.y)
            y_temp[chunk_sum] = y_fit
            if (hasattr(self, '_fit')):
                y_temp[chunk[c]] = self._fit.y[chunk[c]]                
            else:
                y_temp[chunk[c]] = cont[chunk[c]]
            y_try = y_temp[chunk_sum]
            redchi = np.append(redchi,
                               self.redchi(y, dy, y_try, len(y)-fit.nvarys))
            #plt.plot(x, y_try)
        #plt.show()

        return stop

    def fit_auto(self, x=None, line=None, i_max=10, mode=None, **kwargs):
        """ Fit a group of lines 
        
        Given a line, the whole group of adjacient lines is fitted, adding
        components when needed.
        OBSOLETE!
        """

        start = time.time()
        
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

            #print(time.time()-start)
            if i == 1:
                start = voigt_params(self, **kwargs)
            else:
                start = {'z': [], 'N': [], 'b': [], 'btur': []}
            
            # Add associated component
            noneb = dc(self._noneb)
            fit_noneb = noneb.fit_wrap(x, vary, mode, **start)

            #print(time.time()-start)

            # Add generic "nebulium" component
            neb = dc(self._neb)
            fit_neb = neb.fit_wrap(x, vary, mode, **start)

            # Check results and choose the best option
            if ((fit_noneb == None) or (fit_neb == None)):
                stop = True
            else:
                #print(noneb._redchi, neb._redchi)
                if (noneb._redchi <= neb._redchi):
                    self_temp = dc(noneb)
                    fit = fit_noneb
                else:
                    self_temp = dc(neb)
                    fit = fit_neb
                self.__dict__.update(self_temp.__dict__)
                print("(%i) %3.2f;" \
                      % (i, self._redchi), end=" ", flush=True)
                stop = (self._redchi < redchi_thr) \
                       or ((self._redchi<10*redchi_thr) \
                           and (self._aic>aic_old)) \
                       or (i==i_max)
                #or (self._last_add == None) \
                       
                aic_old = self._aic
                redchi_old = self._redchi            
                if (self._redchi < redchi_best): #or 1==1):
                    self_best = dc(self)
                    fit_best = dc(fit)
                    i_best = i
                    redchi_best = self._redchi

                if (stop == False):
                    cont_corr = self.corr_resid(cont_corr)

            #print(time.time()-start)
       
        self = dc(self_best)
        self.__dict__.update(self_best.__dict__)
        fit = fit_best
        print("best chi-squared (%i) %3.2f, %3.2f;" \
              % (i_best, redchi_best, self._aic), end=" ", flush=True)


    def fit_list(self, list_range=None, iter_range=range(5,6), mode=None,
                 plot=True, **kwargs):
        if (list_range is None):
            list_range = range(len(self.t))

        # Read Voigt parameters, if provided
        
        self_temp = dc(self)
        x_arr = self_temp.x
        #i = 0
        group_check = 0
        self._z_list = np.array([])
        self._N_list = np.array([])
        self._b_list = np.array([])
        self._btur_list = np.array([])
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
                    
                    #self.fit_add(x=x_arr[l], i_max=i, mode=mode, **kwargs)
                    self.fit_auto(x=x_arr[l], i_max=i, mode=mode, **kwargs)
                    print("time: %3.2f;" % (time.time()-start), end=" ",
                          flush=True)
                    self._z_list = np.append(self._z_list, self._z_fit) #\
                                   #* self._z_fit.unit
                    self._N_list = np.append(self._N_list, self._N_fit) #\
                                   #* self._N_fit.unit
                    self._b_list = np.append(self._b_list, self._b_fit) #\
                                   #* self._b_fit.unit
                    self._btur_list = np.append(
                        self._btur_list, self._btur_fit) #* self._btur_fit.unit

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

        ### OBSOLETE 
        self._fit_x = self._spec.x[where]
        self._fit_y = self._spec.y[where] / self._cont.y[where]
        self._fit_dy = self._spec.dy[where] / self._cont.y[where]
        ###
        
        self._norm_guess = self.norm(vary=vary)
        #self._norm_guess = self.prof(vary=vary)
        if (prof == 'voigt'):
            if (mode == 'use_old'):
                if (hasattr(self, '_z_fit')):
                    #diff = np.setdiff1d(self._z[self._group[1]], self._z_fit)
                    z_temp = np.append(self._z_fit.value, self._z_add)
                    #print(diff, z_temp)
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
                    #print("fit_prep")
                    #print(z, N, b, btur)
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
            self._prof_guess = self.voigt(z=z, N=N, b=b, btur=btur)

        else:
            raise Exception("Only Voigt profile is supported.")
        self._psf = self.psf()


    def fit_prod(self, fit, prof='voigt'):        
        if (hasattr(self, '_fit') == False):
            self._fit = dc(self._spec)
            self._fit.y = 'nan'
        if (hasattr(self, '_cont') == False):
            self._cont = dc(self._spec)            
            self._cont.y = 'nan'
        if (hasattr(self, '_resid_fit') == False):
            self._resid_fit = dc(self._spec)        
            self._resid_fit.y = 'nan'
        if (hasattr(self, '_resid_cont') == False):
            self._resid_cont = dc(self._spec)            
            self._resid_cont.y = 'nan'
        if (hasattr(self, '_rem') == False):
            self._rem = dc(self._spec)            
            self._rem.y = 'nan'
        if (hasattr(self, '_chunk_sum')):
            where = self._chunk_sum
        else:
            where = np.full(len(self._spec.t), True)
        yunit = self._spec.y.unit
        comp = fit.eval_components(x=self._spec.x[where].value)
        cont = self._cont.y[where]
        slope = cont / np.mean(cont)
        self._fit.y[where] = fit.best_fit * cont
        self._cont.y[where] = comp['cont_'] * cont
        #self._fit.y[where] = fit.best_fit * slope * self._fit.y.unit
        #self._cont.y[where] = comp['cont_'] * slope * self._cont.y.unit
        self._resid_fit.y[where] = self._spec.y[where] - self._fit.y[where]
        self._resid_cont.y[where] = self._spec.y[where] - self._cont.y[where] 
        self._rem.y[where] = self._cont.y[where] + self._resid_fit.y[where]
    #* self._spec.y[where] / self._fit.y[where] #* yunit
        
        if (prof == 'voigt'):
            #print(fit.fit_report())
            #print(fit.params.pretty_print())
            #print(fit.errorbars)
            #print(fit.best_values)
            #print(fit.params)
            #print(fit.params['voigt0_z15138094933394572_btur'].stderr)
            z_tags = [z for z in fit.best_values if z.endswith('_z')]
            N_tags = [N for N in fit.best_values if N.endswith('_N')]
            b_tags = [b for b in fit.best_values if b.endswith('_b')]
            btur_tags = [bt for bt in fit.best_values if bt.endswith('_btur')]

            z_best = np.array([fit.best_values[z] for z in z_tags])
            N_best = np.array([fit.best_values[N] for N in N_tags] )
            b_best = np.array([fit.best_values[b] for b in b_tags]) 
            btur_best = np.array([fit.best_values[bt] for bt \
                                  in btur_tags])

            zerr_best = np.array([fit.params[z].stderr for z in z_tags])
            Nerr_best = np.array([fit.params[N].stderr for N \
                                  in np.sort(N_tags)])
            berr_best = np.array([fit.params[b].stderr for b \
                                  in np.sort(b_tags)])
            bturerr_best = np.array([fit.params[bt].stderr for bt \
                                     in np.sort(btur_tags)])

            z_sort = np.sort(z_best)
            N_sort = N_best[np.argsort(z_best)]
            b_sort = b_best[np.argsort(z_best)]
            btur_sort = btur_best[np.argsort(z_best)]

            zerr_sort = zerr_best[np.argsort(z_best)]
            Nerr_sort = Nerr_best[np.argsort(z_best)]
            berr_sort = berr_best[np.argsort(z_best)]
            bturerr_sort = bturerr_best[np.argsort(z_best)]

            if ('ION' in self.t.colnames):
                sel = np.append(0,
                                np.cumsum([np.size(ion) for ion \
                                in self._t['ION'][self._group[1]]]))[:-1]
            else:
                # Probably obsolete
                sel = range(np.sum(self._group[1]))
            self._z_fit = u.Quantity(z_sort[sel])
            self._N_fit = N_sort[sel] / u.cm**2
            self._b_fit = b_sort[sel] * u.km/u.s
            self._btur_fit = btur_sort[sel] * u.km/u.s

            #print("fit_prod")
            #print(self._z_fit, self._N_fit, self._b_fit, self._btur_fit)
            self._zerr_fit = u.Quantity(zerr_sort[sel])
            self._Nerr_fit = Nerr_sort[sel] / u.cm**2
            self._berr_fit = berr_sort[sel] * u.km/u.s
            self._bturerr_fit = bturerr_sort[sel] * u.km/u.s

            #print(self._Nerr_fit)

            # When new redshift is a duplicate
            """ No action is taken, currently
            if ((hasattr(self, '_last_add')) \
                and (len(self._z_fit) < np.sum(self._group[1]))): 
                #print(self._last_add)
                self._z = np.delete(self._z, self._last_add)
                self._t.remove_row(self._last_add)
                line = self._group[0]
                self.group(line=line)
                #self._group[1] = np.delete(self._group[1], self._last_add)
            #print(len(self._z), len(self._group[1]), len(self._z_fit))
            """
            if (hasattr(self, '_z')):
                self._z[self._group[1]] = self._z_fit

            if ('N' in self._t.colnames):
                self._t['N'][self._group[1]] = self._N_fit
                self._t['B'][self._group[1]] = self._b_fit 
                self._t['BTUR'][self._group[1]] = self._btur_fit
                if (fit.errorbars == True):
                    self._t['DN'][self._group[1]] = self._Nerr_fit
                    self._t['DB'][self._group[1]] = self._berr_fit 
                    self._t['DBTUR'][self._group[1]] = self._bturerr_fit
           
            model = Model(self._spec, syst=self, group=self._group, chunk=self._chunk)
            voigt = model.voigt(self._z_fit, self._N_fit, self._b_fit,
                                self._btur_fit, self._flat.ion)
            self._fit.y[where] = voigt[0].eval(
                voigt[1], x=self._fit.x[where].value) * cont
        else:
            raise Exception("Only Voigt profile is supported.")

        self._redchi = fit.redchi
        self._aic = fit.aic    
        
    def fit_wrap(self, x, vary=False, mode=None, **kwargs):
        """ Model a group of lines an fit them """
        
        group = self.group(x)
        chunk = self.chunk(x)

        self.fit_prep(mode=mode, vary=vary, **kwargs)

        # Create a model
        guess = self.model()

        # Fit the model
        fit = self.fit()

        # Evaluate the fit
        #self.fit_eval(fit)
        
        
        if (hasattr(fit, 'fit_report')):
            self.fit_prod(fit)
        else:
            fit = None
        return fit

    def flatten_z(self):
        """ Create a flattened version of the system, with different entries
        for each ion """

        # This "try" will be removed when fitting methods are moved to "abs"
        try:
            tab = self.t_all
        except:
            tab = self.t
        yunit = tab['Y'].unit
        #yunit = self._linez.dy.unit

        first = True
        for r in tab:
            try:
                # Tuples
                for i in range(len(r['Y'])):
                    if (first == True):
                        #print(r['ION'][i])
                        self._flat = System(x=[r['X']], y=[r['Y'][i]],
                                            xmin=[r['XMIN']], xmax=[r['XMAX']],
                                            dy=[r['DY'][i]], ion=[r['ION'][i]],
                                            N=[r['N']], b=[r['B']],
                                            btur=[r['BTUR']], yunit=yunit)
                        first = False
                    else:
                        row = [r['X'], r['Y'][i], r['XMIN'], r['XMAX'],
                               r['DY'][i], r['ION'][i]]
                        for col in r.colnames[6:]:
                            row.append(r[col])
                        self._flat.t.add_row(row)
            except:
                # Scalars
                self._flat = self            
                            
                            
    def group(self, x=None, line=None, single=False):
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            #line = np.where(abs(self.x.value-x.value) \
            #                == abs(self.x.value-x.value).min())[0][0]
            line = np.where(abs(self.x-x) == abs(self.x-x).min())[0][0]
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
            # This part is needed to add interlopers to the group
            xmin = float('inf') 
            xmax = 0
            # The cycle must run two times because some interloper may affect
            # the other doublet component and would then be missed
            for t in range(2):  
                for r in self._t[sel]:
                    try:
                        # Tuples
                        for i in range(len(r['Y'])):
                            #print(r['ION'][i])
                            xmin = min(xmin,
                                       (1 + r['XMIN']) * dict_wave[r['ION'][i]])
                            xmax = max(xmax,
                                       (1 + r['XMAX']) * dict_wave[r['ION'][i]])
                    except:
                        # Scalars
                        xmin = (1 + r['XMIN']) * dict_wave[r['ION']]
                        xmin = (1 + r['XMAX']) * dict_wave[r['ION']]
                l = 0
                for r in self._t:
                    try:
                        # Tuples 
                        for i in range(len(r['Y'])):
                            x = (1 + r['X']) * dict_wave[r['ION'][i]]
                            if ((x > xmin) and (x < xmax)):
                                sel[l] = True
                    except:
                        # Scalars
                        x = (1 + r['X']) * dict_wave[r['ION']]
                        if ((x > xmin) and (x < xmax)):
                            sel[l] = True
        
                    l += 1
            
        ret = (line, sel)
        self._group = ret
        return ret
        
    def match_z(self, zstart=None, zend=None, ztol=1e-4):
        """ Match redshifts in a list, to define systems """

        if (hasattr(self, '_linez') == False):
            raise Exception("Redshift table must be created before matching.")
        
        
        # Flatten arrays
        # N.B. Y and DY don't need flattening, as they have only one entry
        # per row. They will be repeated to have the shape of the other arrays.
        z = np.ravel(self._linez.x)
        zmin = np.ravel(self._linez.xmin)
        zmax = np.ravel(self._linez.xmax)
        ion = np.ravel(self._linez.ion)
        if (len(self._linez.x.shape) > 1):
            y = np.repeat(self._linez.y, self._linez.x.shape[1])
            dy = np.repeat(self._linez.dy, self._linez.x.shape[1])
        else:
            y = np.ravel(self._linez.y)
            dy = np.ravel(self._linez.dy)

        if (zstart != None and zend != None):
            where = np.logical_and(z > zstart, z < zend)
            z = z[where]
            zmin = zmin[where]
            zmax = zmax[where]            
            ion = ion[where]
            y = y[where]
            dy = dy[where]
            
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
                if (ion_sort[i+1] not in ion_coinc_row):
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
                
        self._z = u.Quantity(np.array(z_coinc))
        syst = System(self.line, self.spec, x=z_coinc, y=y_coinc,
                    xmin=zmin_coinc, xmax=zmax_coinc, dy=dy_coinc,
                    ion=ion_coinc, yunit=y.unit)
        self.__dict__.update(syst.__dict__)

    def norm(self, value=1.0, vary=False):
        """ Normalize continuum """

        model = Model(self._spec, line=self, group=self._group, chunk=self._chunk) 
        norm = model.norm(value, vary)
        if (hasattr(self, '_norm') == False):
            self._norm = dc(self._spec)
        self._norm.y[self._chunk_sum] = norm[0].eval(
            norm[1], x=self._norm.x[self._chunk_sum].value) \
            * self._cont.y[self._chunk_sum] #* self._norm.y[chunk_sum].unit

        return norm 

    def plot(self, group=None, chunk=None, ion=[], figsize=(10,4),
             mode='simple', block=True, **kwargs):
        if (ion == []):
            ion = np.unique(self._flat.ion)
        ion = ion[ion != 'neb']

        n = len(ion)
        try:
            z = self._linez_fit
        except:
            z = self.x

        if (hasattr(self, '_chunk_sum')):
            chunk_sum = self._chunk_sum

        z_neb = np.asarray([z for k in range(len(z)) if z[k] > 30.0])
        # Change this (hardcoded)
        if (len(z_neb) > 0):
            x_neb = (1.0 + z_neb[0]) * dict_wave['neb']
        else:
            x_neb = np.asarray([])
        
        if mode == 'split':
            #figsize = (7,6)
            row = min(n,4)
            col = int(np.ceil(n/4))
            figsize = (col*6, n*3.5)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            grid = gs(row,col)
            for p in range(n):
                zmin = float('inf')
                zmax = 0
                if (group is not None):
                    t = Table(self._t[group[1]])
                    for l in range(len(t)):
                        if (ion[p] in t['ION'][l]):
                            zmin = min(zmin, t['XMIN'][l])
                            zmax = max(zmax, t['XMAX'][l])
                        #zmin = min(zmin, t['XMIN'][l])
                        #zmax = max(zmax, t['XMAX'][l])
                else:
                    zmin = np.min(self.xmin)
                    zmax = np.max(self.xmax)
                ax = fig.add_subplot(grid[p%4, int(np.floor(p/4))])
                ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
                """
                if (hasattr(self._spec, '_orig')):
                    spec = dc(self._spec._orig)
                else:
                    spec = dc(self._spec)
                """
                spec = dc(self._spec)
                spec.to_z([ion[p]])
                if self._line != None:
                    line = dc(self._line)
                    line.to_z([ion[p]])
                else:
                    line = self
                ax.plot(spec.x, spec.y, c='black', lw=1.0)
                #ax.plot(spec.x, spec.dy, c='r', lw=1.0)
                #ax.plot(spec.x, -spec.dy, c='r', lw=1.0)
                if (chunk is not None):
                    if (hasattr(self, '_norm')):
                        norm = dc(self._norm)
                        norm.to_z([ion[p]])
                        """
                        for c in range(1, len(chunk)):
                            ax.plot(norm.x[chunk[c]], norm.y[chunk[c]], c='y',
                                    lw=1.0, linestyle=':')
                        """
                        ax.plot(norm.x[chunk_sum], norm.y[chunk_sum], c='y',
                                lw=1.0, linestyle=':')
                    if (hasattr(self, '_voigt')):
                        voigt = dc(self._voigt)
                        voigt.to_z([ion[p]])
                        """
                        for c in range(1, len(chunk)):
                            ax.plot(voigt.x[chunk[c]], voigt.y[chunk[c]], c='g',
                                    lw=1.0, linestyle=':')
                        """
                        ax.plot(voigt.x[chunk_sum], voigt.y[chunk_sum], c='g',
                                lw=1.0, linestyle=':')
                    if (hasattr(self, '_cont')):
                        cont = dc(self._cont)
                        cont.to_z([ion[p]])
                        """
                        for c in range(1, len(chunk)):
                            ax.plot(cont.x[chunk[c]], cont.y[chunk[c]], c='y')
                        """
                        ax.plot(cont.x[chunk_sum], cont.y[chunk_sum], c='y')
                    if (hasattr(self, '_fit')):
                        fit = dc(self._fit)
                        fit.to_z([ion[p]])
                        """
                        for c in range(1, len(chunk)):
                            ax.plot(fit.x[chunk[c]], fit.y[chunk[c]], c='g')
                        """
                        ax.plot(fit.x[chunk_sum], fit.y[chunk_sum], c='g')
                    """    
                    if (hasattr(self, '_resid_fit')):
                        resid_fit = dc(self._resid_fit)
                        resid_fit.to_z([ion[p]])
                        ax.plot(resid_fit.x[chunk_sum],
                                resid_fit.y[chunk_sum], c='b', lw=1.0)
                    """    
                    """
                    if (hasattr(self, '_rem')):
                        rem = dc(self._rem)
                        rem.to_z([ion[p]])
                        ax.plot(rem.x[chunk_sum], rem.y[chunk_sum], c='b',
                                lw=1.0)
                    """
                #ax.scatter(line.x, line.y, c='b')
                for comp in z:
                    ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='black')
                if (len(x_neb) > 0):
                    for comp_neb in x_neb/dict_wave[ion[p]] - 1.0:
                        ax.axvline(x=comp_neb, ymin=0.65, ymax=0.85, color='r')
                if ((p+1) % row != 0):
                    pass
                    #ax.set_xticks([], [])
                else:
                    ax.set_xlabel("Redshift")
                ax.set_xlim(zmin, zmax)
                ax.set_ylim(np.min(spec.y[chunk_sum].value),
                            np.max([np.max(cont.y[chunk_sum].value),
                                    np.max(spec.y[chunk_sum].value)]))
                ax.text(0.5, 0.92, ion[p], horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes,
                        fontsize=12)
            if (hasattr(self, '_fit')):
                fig.suptitle("Reduced chi-squared: %3.2f" % (self._redchi),
                             fontsize=10)
        elif mode == 'compare':
            if (chunk is not None):
                zmin = np.min(self.xmin[group[1]])
                zmax = np.max(self.xmax[group[1]])
            else:
                zmin = np.min(self.xmin)
                zmax = np.max(self.xmax)
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
                ax.plot(spec.x, spec.y, lw=1.0)
                ax.scatter(line.x, line.y, marker='+')
            text = ', '.join(str(p) for p in ion)
            for comp in z:
                ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='black', lw=3.0,
                           linestyle='--')

            ax.text(0.5, 0.92, text, horizontalalignment="center",
                    verticalalignment="center", transform=ax.transAxes,
                    fontsize=12)
        elif mode == 'simple':
            zmin = np.min(self.xmin)
            zmax = np.max(self.xmax)
            grid = gs(1,1)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            ax = fig.add_subplot(grid[:,:])
            ax.set_xlabel("Redshift")
            ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
            spec = dc(self._spec)
            line = dc(self._line)
            spec.to_z([ion[0]])
            line.to_z([ion[0]])
            ax.set_xlim(zmin, zmax)
            ax.plot(spec.x, spec.y, c='black', lw=1.0)
            ax.plot(spec.x, spec.dy, c='r', lw=1.0)
            if (chunk is not None):
                if (hasattr(self, '_norm')):
                    norm = dc(self._norm)
                    norm.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(norm.x[chunk[c]], norm.y[chunk[c]], c='y',
                                lw=1.0, linestyle=':')
                if (hasattr(self, '_voigt')):
                    voigt = dc(self._voigt)
                    voigt.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(voigt.x[chunk[c]], voigt.y[chunk[c]], c='g',
                                lw=1.0, linestyle=':')
                if (hasattr(self, '_cont')):
                    cont = dc(self._cont)
                    cont.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(cont.x[chunk[c]], cont.y[chunk[c]], c='y')
                if (hasattr(self, '_fit')):
                    fit = dc(self._fit)
                    fit.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(fit.x[chunk[c]], fit.y[chunk[c]], c='g')
                if (hasattr(self, '_rem')):
                    rem = dc(self._rem)
                    rem.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(rem.x[chunk[c]], rem.y[chunk[c]], c='b', lw=1.0)
            ax.scatter(line.x, line.y, marker='+')
            text = ion[0] + ' redshifts'

            ax.text(0.5, 0.92, text, horizontalalignment="center",
                    verticalalignment="center", transform=ax.transAxes,
                    fontsize=12)
           
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
        grid.update(hspace=0.2)
        if block is True:
            plt.show()

    def psf(self):
        """ Model the instrumental PSF """
        
        model = Model(self._spec, line=self, group=self._group,
                      chunk=self._chunk)
        psf = model.psf()

        return psf   

    def psf2(self, resol):
        """ Model the instrumental PSF """
        
        resol_arr = np.ones(np.sum(self._chunk_sum)) * resol
        model = Model(self._spec, line=self, group=self._group,
                      chunk=self._chunk)
        psf = model.psf2(resol_arr)
        return psf   
    
  #  def redchi(self, model_param, nvarys):
    def redchi(self, y, dy, y_fit, dof):
        """
        model = model_param[0]
        param = model_param[1]
        x = self._spec.x[self._chunk_sum]
        y = self._spec.y[self._chunk_sum]        
        dy = self._spec.dy[self._chunk_sum]
        ndata = len(x)
        mod = model.eval(param, x=x.value)
        ret = np.sum(((mod-y.value)/dy.value)**2) / (ndata-nvarys)

        """
        ret = np.sum(((y_fit-y)/dy)**2) / dof
        return ret

    def save(self, name):
        """
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', array=self._spec.xmin),
             fits.Column(name='XMAX', format='E', array=self._spec.xmax),
             fits.Column(name='X', format='E', array=self._spec.x),
             fits.Column(name='Y', format='E', array=self._spec.y),
             fits.Column(name='Y_FIT', format='E', array=self._fit.y),
             fits.Column(name='Y_REM', format='E', array=self._rem.y),
             fits.Column(name='DY', format='E', array=self._spec.dy),
             fits.Column(name='GROUP', format='I', array=self._spec.group),
             fits.Column(name='RESOL', format='E', array=self._spec.resol)]) 
        hdu.writeto(name + '_syst_spec.fits', overwrite=True)
        """
        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='X', format='E', array=self.x),
             fits.Column(name='XMIN', format='E', array=self.xmin),
             fits.Column(name='XMAX', format='E', array=self.xmax),
             #fits.Column(name='Y', format='E', array=self._t['Y']),
             #fits.Column(name='DY', format='E', array=self._t['DY']),
             #fits.Column(name='ION', format='I', array=self.ion),
             fits.Column(name='N', format='E', array=self._t['N']),
             fits.Column(name='B', format='E', array=self._t['B']),
             fits.Column(name='BTUR', format='E', array=self._t['BTUR'])]) 
        hdu.writeto(name + '_syst.fits', overwrite=True)

        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', array=self._fit.xmin),
             fits.Column(name='XMAX', format='E', array=self._fit.xmax),
             fits.Column(name='X', format='E', array=self._fit.x),
             fits.Column(name='Y', format='E', array=self._fit.y),
             fits.Column(name='DY', format='E', array=self._fit.dy),
             fits.Column(name='GROUP', format='I', array=self._fit.group),
             fits.Column(name='RESOL', format='E', array=self._fit.resol)])
        hdu.writeto(name + '_syst_fit.fits', overwrite=True)

        hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='XMIN', format='E', array=self._rem.xmin),
             fits.Column(name='XMAX', format='E', array=self._rem.xmax),
             fits.Column(name='X', format='E', array=self._rem.x),
             fits.Column(name='Y', format='E', array=self._rem.y),
             fits.Column(name='DY', format='E', array=self._rem.dy),
             fits.Column(name='GROUP', format='I', array=self._rem.group),
             fits.Column(name='RESOL', format='E', array=self._rem.resol)])
        hdu.writeto(name + '_syst_rem.fits', overwrite=True)

        
    def unabs(self):
        """ Remove lines """

        model = Model(self._spec, syst=self, group=self._group, chunk=self._chunk) 
        unabs = model.unabs()
        if (hasattr(self, '_unabs') == False):
            self._unabs = dc(self._spec)
        for c in range(1, len(self._chunk)):
            self._unabs.y[self._chunk[c]] = unabs[2*(c-1)].eval(
                unabs[2*c-1], x=self._unabs.x[self._chunk[c]].value) \
                * self._unabs.y[self._chunk[c]].unit
    
        return unabs

    def voigt(self, z=[], N=[], b=[], btur=[]):

        sumlen = len(z) + len(N) + len(b) + len(btur)
        if ((z != []) and (sumlen % len(z) != 0)):
            raise Exception("Parameter arrays must have the same length.")

        model = Model(self._spec, syst=self, group=self._group, chunk=self._chunk)
        if (z == []):
            z = self._z[self._group[1]]
            if (hasattr(self, '_norm')):
                N = model.N_guess(self._norm, ion=self._flat.ion)
            else:
                N = model.N_guess(self._unabs, ion=self._flat.ion)
            if (hasattr(self, '_unabs')):
                cont = self._unabs
            elif (hasattr(self, '_norm')):
                cont = self._norm
            else:
                raise Exception("Continuum not found.")
            N = model.N_guess(cont, ion=self._flat.ion)
            b = np.full(len(self.x[self._group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(self.x[self._group[1]]), voigt_def['btur']) \
                   * u.km / u.s
        else:
            for val in N:
                if (val == voigt_def['N']):
                    N = model.N_guess(self._norm, ion=self._flat.ion)
            
        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)

        ion = np.unique(self._flat.ion)
        #print("voigt")
        #print(z, N, b, btur)
        voigt = model.voigt(z, N, b, btur, ion)
        """
        for c in range(1, len(chunk)):
            if (c == 1):
                chunk_sum = dc(chunk[c])
            else: 
                chunk_sum += chunk[c]
        """
        #print("voigt")
        #voigt[1].pretty_print()
        self._voigt.y[self._chunk_sum] = voigt[0].eval(voigt[1], 
            x=self._voigt.x[self._chunk_sum].value) * self._voigt.y.unit
        if (hasattr(self, '_norm')):
            self._voigt.y[self._chunk_sum] = self._voigt.y[self._chunk_sum] \
                                       * self._norm.y[self._chunk_sum].value
        else:
            self._voigt.y[self._chunk_sum] = self._voigt.y[self._chunk_sum] \
                                       * self._unabs.y[self._chunk_sum].value    
            
        self._z_arr = dc(model._z)
        self._N_arr = dc(model._N)
        self._b_arr = dc(model._b)
        self._btur_arr = dc(model._btur)

        return voigt
