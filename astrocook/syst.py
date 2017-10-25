from . import Line, Model
from .utils import convolve, convolve2, dict_doubl, dict_wave, voigt_def
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
                 doubl=[],
                 N=[],
                 b=[],
                 btur=[],
                 yunit=None,
                 meta=None,
                 dtype=float):
        """ Constructor for the Syst class """ 
        
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
            ion = np.stack((np.full(len_ion, dict_doubl[doubl][0]),
                            np.full(len_ion, dict_doubl[doubl][1]))).T
                           
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
        if isinstance(value, Syst):
            self._linez = value
        else:
            raise Exception("Redshift list has a wrong format.")


# Methods

    #def add_assoc(self, ion, z=None, zmin, zmax):
    #    self.t.add_row([z, y, zmin, zmax, dy, ion])

    def corr_resid(self, group, cont_corr):
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
            tab = Table(self.t[group[1]][ion_arr[p] in 'ION'])
            z_cen[p] = tab['X']

        where = (abs(z_arr-z_cen) == abs(z_arr-z_cen).min())
        z_ion = z_arr[where][0]

        ion_where = (abs(z_ion-self.x) == abs(z_ion-self.x).min())
        ion = self.ion[ion_where][0]
        zmin_ion = self.xmin[ion_where][0] 
        zmax_ion = self.xmax[ion_where][0]
        #y_ion = np.empty(len(ion))
        #dy_ion = np.empty(len(ion))        
        for i in range(len(ion)):
            spec = dc(self._spec)
            spec.to_z([ion[i]])
            #y_ion[i] = np.interp(z_ion, spec.x, spec.y)
            #dy_ion[i] = np.interp(z_ion, spec.x, spec.dy)
            if (i == 0):
                y_ion = (np.interp(z_ion, spec.x, spec.y),)
                dy_ion = (np.interp(z_ion, spec.x, spec.dy),)
            else:
                y_ion += (np.interp(z_ion, spec.x, spec.y),)
                dy_ion += (np.interp(z_ion, spec.x, spec.dy),)
                

        
        z_neb = x / dict_wave[neb] - 1.0
        line = dc(self._line)
        line.to_z([neb])
        neb_where = abs(z_neb-line.x) == abs(z_neb-line.x).min()
        zmin_neb = line.xmin[neb_where][0] 
        zmax_neb = line.xmax[neb_where][0]
        spec = dc(self._spec)
        spec.to_z([neb])
        y_neb = np.interp(x.value, self._spec.x, self._spec.y) \
                * self._spec.yunit
        dy_neb = np.interp(x.value, self._spec.x, self._spec.dy) \
                 * self._spec.yunit

        
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

        syst = Syst(self.line, self.spec, x=z_deflat, y=y_deflat,
                    xmin=zmin_deflat, xmax=zmax_deflat, dy=dy_deflat,
                    ion=ion_deflat, N=N_deflat, b=b_deflat, btur=btur_deflat,
                    yunit=yunit)
        self.__dict__.update(syst.__dict__)

    def find(self, ztol=1e-4):
        self.create_z()
        self.match_z(ztol)
        self.flatten_z()

        
    def flatten_z(self):
        """ Create a flattened version of the system, with different entries
        for each ion """

        yunit = self._linez.dy.unit
        first = True
        for r in self.t:
            for i in range(len(r['Y'])):
                if (first == True):
                    self._flat = Syst(x=[r['X']], y=[r['Y'][i]],
                                      xmin=[r['XMIN']], xmax=[r['XMAX']],
                                      dy=[r['DY'][i]], ion=[r['ION'][i]],
                                      N=[r['N']], b=[r['B']], btur=[r['BTUR']],
                                      yunit=yunit)
                    first = False
                else:
                    row = [r['X'], r['Y'][i], r['XMIN'], r['XMAX'], r['DY'][i],
                           r['ION'][i]]
                    for col in r.colnames[6:]:
                        row.append(r[col])
                    self._flat.t.add_row(row)
                    
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
                    for i in range(len(r['Y'])):
                        xmin = min(xmin,
                                   (1 + r['XMIN']) * dict_wave[r['ION'][i]])
                        xmax = max(xmax,
                                   (1 + r['XMAX']) * dict_wave[r['ION'][i]])
                l = 0
                for r in self._t:
                    for i in range(len(r['Y'])):
                        x = (1 + r['X']) * dict_wave[r['ION'][i]]
                        if ((x > xmin) and (x < xmax)):
                            sel[l] = True
                    l += 1
            
        ret = (line, sel)
        self._group = ret
        return ret
        
    def match_z(self, ztol=1e-4):
        """ Match redshifts in a list, to define systems """

        if (hasattr(self, '_linez') == False):
            raise Exception("Redshift table must be created before matching.")
        
        
        # Flatten arrays
        # N.B. Y and DY don't need flattening, as they have only one entry
        # per row. They will be repeated to have the shape of the other arrays.
        z = np.ravel(self._linez.x)
        y = np.repeat(self._linez.y, self._linez.x.shape[1])        
        zmin = np.ravel(self._linez.xmin)
        zmax = np.ravel(self._linez.xmax)
        dy = np.repeat(self._linez.dy, self._linez.x.shape[1])
        ion = np.ravel(self._linez.ion)

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
        syst = Syst(self.line, self.spec, x=z_coinc, y=y_coinc,
                    xmin=zmin_coinc, xmax=zmax_coinc, dy=dy_coinc,
                    ion=ion_coinc, yunit=y.unit)
        self.__dict__.update(syst.__dict__)

    def norm(self, group, chunk, value=1.0, vary=False):
        """ Normalize continuum """

        model = Model(self._spec, line=self, group=group, chunk=chunk) 
        norm = model.norm(value, vary)
        if (hasattr(self, '_norm') == False):
            self._norm = dc(self._spec)
        for c in range(1, len(chunk)):
            #self._norm.y[chunk[c]] = norm[2*(c-1)].eval(
            #    norm[2*c-1], x=self._norm.x[chunk[c]].value) \
            #    * self._cont.y[chunk[c]] * self._norm.y[chunk[c]].unit 
            if (c == 1):
                chunk_sum = dc(chunk[c])
            else: 
                chunk_sum += chunk[c]
        self._norm.y[chunk_sum] = norm[0].eval(
            norm[1], x=self._norm.x[chunk_sum].value) \
            * self._cont.y[chunk_sum] #* self._norm.y[chunk_sum].unit 

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
                line = dc(self._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
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
                ax.scatter(line.x, line.y, c='b')
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
                            np.max(spec.y[chunk_sum].value))
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
            text = ion[0] + 'redshifts'

            ax.text(0.5, 0.92, text, horizontalalignment="center",
                    verticalalignment="center", transform=ax.transAxes,
                    fontsize=12)
           
        grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.97])
        grid.update(hspace=0.2)
        if block is True:
            plt.show()

    def psf(self, group, chunk, resol):
        """ Model the instrumental PSF """
        
        model = Model(self._spec, line=self, group=group, chunk=chunk)
        psf = model.psf(resol)

        return psf   

    def psf2(self, group, chunk, resol):
        """ Model the instrumental PSF """
        
        for c in range(1, len(chunk)):
            if (c == 1):
                chunk_sum = dc(chunk[c])
            else:
                chunk_sum += chunk[c]
        resol_arr = np.ones(np.sum(chunk_sum)) * resol
        model = Model(self._spec, line=self, group=group, chunk=chunk)
        psf = model.psf2(resol_arr)
        return psf   
    
    def redchi(self, model_param, nvarys):
        model = model_param[0]
        param = model_param[1]
        for c in range(1, len(self._chunk)):
            if (c == 1):
                chunk_sum = dc(self._chunk[c])
            else:
                chunk_sum += self._chunk[c]
        x = self._spec.x[chunk_sum]
        y = self._spec.y[chunk_sum]        
        dy = self._spec.dy[chunk_sum]
        ndata = len(x)
        mod = model.eval(param, x=x.value)
        ret = np.sum(((mod-y.value)/dy.value)**2) / (ndata-nvarys)
        return ret

    def save(self, name):
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

    def voigt(self, group=None, chunk=None, z=[], N=[], b=[], btur=[]):

        sumlen = len(z) + len(N) + len(b) + len(btur)
        if ((z != []) and (sumlen % len(z) != 0)):
            raise Exception("Parameter arrays must have the same length.")

        model = Model(self._spec, syst=self, group=group, chunk=chunk)
        if (z == []):
            #z = self.x[group[1]]
            z = self._z[group[1]]
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
            b = np.full(len(self.x[group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(self.x[group[1]]), voigt_def['btur']) \
                   * u.km / u.s

        else:
            for val in N:
                if (val == voigt_def['N']):
                    N = model.N_guess(self._norm, ion=self._flat.ion)
            
        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)

        ion = np.unique(self._flat.ion)
        voigt = model.voigt(z, N, b, btur, ion)
        for c in range(1, len(chunk)):
            """
            self._voigt.y[chunk[c]] = voigt[2*(c-1)].eval(
                voigt[2*c-1], x=self._voigt.x[chunk[c]].value) \
                * self._voigt.y.unit
            if (hasattr(self, '_norm')):
                self._voigt.y[chunk[c]] = self._voigt.y[chunk[c]] \
                                          * self._norm.y[chunk[c]].value
            else:
                self._voigt.y[chunk[c]] = self._voigt.y[chunk[c]] \
                                          * self._unabs.y[chunk[c]].value
            """
            if (c == 1):
                chunk_sum = dc(chunk[c])
            else: 
                chunk_sum += chunk[c]
        self._voigt.y[chunk_sum] = voigt[0].eval(voigt[1], 
            x=self._voigt.x[chunk_sum].value) * self._voigt.y.unit
        if (hasattr(self, '_norm')):
            self._voigt.y[chunk_sum] = self._voigt.y[chunk_sum] \
                                       * self._norm.y[chunk_sum].value
        else:
            self._voigt.y[chunk_sum] = self._voigt.y[chunk_sum] \
                                       * self._unabs.y[chunk_sum].value    

        self._z_arr = dc(model._z)
        self._N_arr = dc(model._N)
        self._b_arr = dc(model._b)
        self._btur_arr = dc(model._btur)

        return voigt
