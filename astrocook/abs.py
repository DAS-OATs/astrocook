from . import Line, Model, Syst
from .utils import convolve, convolve2, dict_doubl, dict_wave, redchi_thr, voigt_def
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

class Abs(Syst):

    def __init__(self,
                 source,
                 ion=[],
                 doubl=[],
                 N=[],
                 b=[],
                 btur=[],
                 meta=None):
        """ Constructor for the Abs class """ 
        
        """
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
        """

        # Spectrum
        self._spec = dc(source._spec)

        # Line list
        # If source is a system
        if (hasattr(source, 'line')):  
            self._line = dc(source._line) 
            self._syst = dc(source)
            source.t.remove_columns(['N', 'B', 'BTUR', 'DN', 'DB', 'DBTUR'])
        # If source is a line list
        else:                          
            self._line = dc(source)
            self._syst = None
        self._source = dc(source) 


        
        if (hasattr(source, '_cont')):
            self._precont = dc(source._precont)
            self._cont = dc(source._cont)
        if (hasattr(source, '_minima')):
            self._minima = dc(source._minima)
        if (hasattr(source, '_maxima')):
            self._maxima = dc(source._maxima)            


        # Ion list
        """
        if (line is not None):
            len_ion = len(line._t)
        if (x != []):
            len_ion = len(x)
        """
        len_ion = len(source.t)
        if ((ion == []) and (doubl != [])):
            ion = np.stack((np.full(len_ion, dict_doubl[doubl][0]),
                            np.full(len_ion, dict_doubl[doubl][1]))).T
                           
        if ((ion == []) or (ion == 'Ly_a')):
            ion = np.full(len_ion, 'Ly_a')           

        if (N == []):
            N = np.full(len_ion, float('nan')) 
            b = np.full(len_ion, float('nan')) 
            btur = np.full(len_ion, float('nan'))
        dN = np.full(len_ion, float('nan')) 
        db = np.full(len_ion, float('nan')) 
        dbtur = np.full(len_ion, float('nan'))
        
        if (hasattr(source, '_ion')):
            self._ion = source._ion
        else:
            self._ion = ion

        # System list
        data = ()
        data_all = ()

        yunit = source.yunit
        col_x = Column(np.asarray(dc(source.t['X'])), name='X')
        col_y = Column(np.asarray(dc(source.t['Y'])), name='Y', unit=yunit)
        col_xmin = Column(np.asarray(dc(source.t['XMIN'])), name='XMIN')
        col_xmax = Column(np.asarray(dc(source.t['XMAX'])), name='XMAX')
        col_dy = Column(np.asarray(dc(source.t['DY'])), name='DY', unit=yunit)
        data_all = (col_x, col_y, col_xmin, col_xmax, col_dy)
            
        #if (ion != []):
        col_ion = Column(np.asarray(dc(self._ion)), name='ION')
        col_N = Column(np.asarray(dc(N)), name='N', unit=1/u.cm**2)
        col_b = Column(np.asarray(dc(b)), name='B', unit=u.km/u.s)
        col_btur = Column(np.asarray(dc(btur)), name='BTUR',
                          unit=u.km/u.s)
        col_dN = Column(np.asarray(dc(dN)), name='DN', unit=1/u.cm**2)
        col_db = Column(np.asarray(dc(db)), name='DB', unit=u.km/u.s)
        col_dbtur = Column(np.asarray(dc(dbtur)), name='DBTUR',
                           unit=u.km/u.s)
        data = (col_ion, col_N, col_b, col_btur, col_dN, col_db, col_dbtur)
        data_all += data
        
        if data is ():
            data = None
            
        if (meta is None):
            meta = {}

        self._t = Table(data=data, masked=True, meta=meta)
        self._t_all = Table(data=data_all, masked=True, meta=meta)
        """
        if (y != []):
            self._t['Y'].unit = yunit
        if (dy != []):
            self._t['DY'].unit = yunit
        """    
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

    @property
    def t(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t

    @property
    def t_all(self):
        if self._use_good:
            return self._t_all[self._igood]
        else:
            return self._t_all

# Methods

    #def add_assoc(self, ion, z=None, zmin, zmax):
    #    self.t.add_row([z, y, zmin, zmax, dy, ion])

    #"""
    def add(self, group, cont_corr, ion=None):

        source = self
        #print(source.t)
        neb = 'neb'

        
        where = source._chunk_sum
        resid_norm = np.full(len(source._resid_fit.y.value), 
                             np.max(source._resid_fit.y[where]\
                                    /source._resid_fit.dy[where]) * 10)
        resid_norm[where] = source._resid_fit.y[where]/source._resid_fit.dy[where]
        x = source._resid_fit.x[where][np.argmin(resid_norm[where])]
        y = np.interp(x.value, self._spec.x, self._spec.y) * self._spec.yunit
        dy = np.interp(x.value, self._spec.x, self._spec.dy) * self._spec.yunit

        ion_arr = np.unique(source._flat.ion)
        n = len(ion_arr)
        z_arr = np.empty(n)        
        z_cen = np.empty(n)        

        #xmin_ion = float('inf') * u.nm
        #xmax_ion = 0 * u.nm
        
        for p in range(n):
            z_arr[p] = x / dict_wave[ion_arr[p]] - 1.0
            tab = Table(source.t_all[group[1]][ion_arr[p] in 'ION'])
            z_cen[p] = tab['X']

        where = (abs(z_arr-z_cen) == abs(z_arr-z_cen).min())
        z_ion = z_arr[where][0]

        ion_where = (abs(z_ion-source.t_all['X']) == abs(z_ion-source.t_all['X']).min())
        ion = source.ion[ion_where][0]
        zmin_ion = source.t_all['XMIN'][ion_where][0] 
        zmax_ion = source.t_all['XMAX'][ion_where][0]
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
        self._noneb.t_all.add_row([z_ion, y_ion, zmin_ion, zmax_ion, dy_ion, 
                               ion, float('nan'), float('nan'), float('nan'),
                               float('nan'), float('nan'), float('nan')])
        self._noneb.t_all.sort('X')  # This gives an annoying warning
        self._noneb._z = np.unique(np.append(source._z.value, z_ion)) 
        self._noneb._z = np.append(source._z.value, z_ion) 
        self._noneb._z.sort()
        self._noneb._z *= source._z.unit
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
        self._neb._t_all = dc(self._neb.t)
        self._neb.t_all.sort('X')  # This gives an annoying warning
        self._neb._z = np.unique(np.append(source._z.value, z_neb.value))
        #self._neb._z = np.append(source._z.value, z_neb.value)
        self._neb._z.sort()
        self._neb._z *= source._z.unit
        self._neb._last_add = np.where(self._neb._z == z_neb)[0][0]
        #else:
        #    self._neb._last_add = None

        self._noneb._source._t = self._noneb.t_all['X', 'Y', 'XMIN', 'XMAX',
                                                  'DY', 'ION']
        self._noneb._t = self._noneb.t_all['ION', 'N', 'B', 'BTUR', 'DN', 'DB',
                                          'DBTUR']
        self._neb._source._t = self._neb.t_all['X', 'Y', 'XMIN', 'XMAX',
                                               'DY', 'ION']
        self._neb._t = self._neb.t_all['ION', 'N', 'B', 'BTUR', 'DN', 'DB',
                                       'DBTUR']
        #print(self._source.t)
        """
        print(self._noneb.t_all)
        print(self._noneb._source.t)
        print(self._noneb.t)
        print(self._neb.t_all)
        print(self._neb._source.t)
        print(self._neb.t)
        """
        return cont_corr

    def chunk(self, x=None, line=None, single=False):  # Chunk must be shifted to the system z
        source = self._source
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            #line = np.where(abs(source.x-x.value) \
            #                == abs(source.x-x.value).min())[0][0]
            line = np.where(abs(source.x-x) == abs(source.x-x).min())[0][0]
        if ((x is None) and (line >= len(source._t))):
            raise Exception("Line number is too large.")

        try:  # When ION has different sizes in different rows
            ion = np.unique(np.asarray(np.sum(source.ion)))
        except:
            ion = np.unique(source.ion)
        n = len(ion)
        iter = range(len(source._t))
        ret = (line,)
        for p in range(n):
            sel = self._spec.t['X'] < 0.0
            spec = dc(self._spec)
            spec.to_z([ion[p]])
            for row in source.t[self.group(line=line, single=single)[1]]:
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


    def fit_add(self, x=None, line=None, i_max=10, mode=None):
        """ Fit a line group, automatically adding components """
        source = self._source
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (line is not None):
            if (x is not None):
                warnings.warn("x will be used; line will be disregarded.")
            x = source.x[line]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        stop = False
        aic_old = float('inf')
        redchi_old = float('inf')        
        redchi_best = float('inf')
        i = 0
        i_best = 1
        cont_corr = 1.0
        vary = True #False
        self._noneb = dc(source)
        self._neb = dc(source)
        self._last_add = 0.0
        while (stop == False):
            i += 1

            # "Non-nebulium" component
            noneb = dc(self._noneb)
            fit_noneb = noneb.fit_wrap(x, vary, mode)

            # Generic "nebulium" component
            neb = dc(self._neb)
            fit_neb = neb.fit_wrap(x, vary, mode)
            
            if ((fit_noneb == None) or (fit_neb == None)):
                stop = True
            else:
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
                    cont_corr = self.add(self._group, cont_corr)

        
        self = dc(self_best)
        self.__dict__.update(self_best.__dict__)
        fit = fit_best
        print("best chi-squared (%i) %3.2f, %3.2f;" \
              % (i_best, redchi_best, self._aic), end=" ", flush=True)

    def fit_all(self, list_range=None, iter_range=range(5,6), mode=None,
                 plot=True):
        if (list_range is None):
            list_range = range(len(self.t))

        self_temp = dc(self)
        x_arr = self._source.x
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
                    self.fit_add(x=x_arr[l], i_max=i, mode=mode)
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
                        self._source._chunk_sum = self._chunk_sum
                        self._source._fit = self._fit
                        self._source._redchi = self._redchi
                        self._source.plot(self._group, self._chunk, mode='split')
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

            #print(self.t)
            #print(self._flat.t)
            self.flatten_z()
            self._prof_guess = self.voigt(self._group, self._chunk, z=z, N=N,
                                          b=b, btur=btur)

        else:
            raise Exception("Only Voigt profile is supported.")
        self._psf = self.psf(self._group, self._chunk, self._source._resol)

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
            N_best = np.array([fit.best_values[N] for N in np.sort(N_tags)] )
            b_best = np.array([fit.best_values[b] for b in np.sort(b_tags)]) 
            btur_best = np.array([fit.best_values[bt] for bt \
                                  in np.sort(btur_tags)])

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
                                np.cumsum([len(ion) for ion \
                                in self._t['ION'][self._group[1]]]))[:-1]
            else:
                sel = range(np.sum(self._group[1]))
            
            self._z_fit = u.Quantity(z_sort[sel])
            self._N_fit = N_sort[sel] / u.cm**2
            self._b_fit = b_sort[sel] * u.km/u.s
            self._btur_fit = btur_sort[sel] * u.km/u.s

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
           
        else:
            raise Exception("Only Voigt profile is supported.")

        self._redchi = fit.redchi
        self._aic = fit.aic    

    def fit_wrap(self, x, vary=False, mode=None):
        source = self #._source
        group = source.group(x)
        chunk = source.chunk(x)
        source.fit_prep(mode=mode, vary=vary)
        guess = source.model()
        #print(source._cont.y)
        fit = source.fit()
        if (hasattr(fit, 'fit_report')):
            source.fit_prod(fit)
        else:
            fit = None
        self._redchi = source._redchi
        self._aic = source._aic
        return fit

    def group(self, x=None, line=None, single=False):
        source = self._source
        
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            #line = np.where(abs(self.x.value-x.value) \
            #                == abs(self.x.value-x.value).min())[0][0]
            line = np.where(abs(source.x-x) ==
                            abs(source.x-x).min())[0][0]
        if ((x is None) and (line >= len(source._t))):
            raise Exception("line is too large.")

        if (single == True):
            sel = np.full(len(source._t), False)
            sel[line] = True
        else:
            iter = range(len(source._t))
            source._t.sort('X')  # This gives a warning on a masked table
            #source._t.group_by('X')
            groups = np.append([0], np.cumsum(source._t['XMAX'][:-1] <
                                              source._t['XMIN'][1:])) 
            sel = np.array([groups[l] == groups[line] for l in iter])
            
            # This part is needed to add interlopers to the group
            xmin = float('inf') 
            xmax = 0
            # The cycle must run two times because some interloper may affect
            # the other doublet component and would then be missed
            for t in range(2):
                for (r, s) in zip(source._t[sel], self._t[sel]):
                    for i in range(np.size(r['Y'])):
                        if (np.size(s['ION']) == 1):
                            try:
                                ion = np.array(s['ION'])[0]
                            except:
                                ion = s['ION']
                        else:
                            ion = s['ION'][i]
                        xmin = min(xmin,
                                   (1 + r['XMIN']) * dict_wave[ion])
                        xmax = max(xmax,
                                   (1 + r['XMAX']) * dict_wave[ion])
                l = 0
                for (r, s) in zip(source._t, self._t):
                    for i in range(np.size(r['Y'])):
                        if (np.size(s['ION']) == 1):
                            try:
                                ion = np.array(s['ION'])[0]
                            except:
                                ion = s['ION']
                        else:
                            ion = s['ION'][i]
                        x = (1 + r['X']) * dict_wave[ion]
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

    """
    def plot(self, group=None, chunk=None, ion=[], figsize=(10,4),
             mode='simple', block=True, **kwargs):

        source = self._source
        if (ion == []):
            ion = np.unique(source._flat.ion)
        ion = ion[ion != 'neb'] 
        n = len(ion)
        try:
            z = source._linez_fit
        except:
            z = source.x

        if (hasattr(source, '_chunk_sum')):
            chunk_sum = source._chunk_sum

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
                    t = Table(source._t[group[1]])
                    for l in range(len(t)):
                        if (ion[p] in t['ION'][l]):
                            zmin = min(zmin, t['XMIN'][l])
                            zmax = max(zmax, t['XMAX'][l])
                        #zmin = min(zmin, t['XMIN'][l])
                        #zmax = max(zmax, t['XMAX'][l])
                else:
                    zmin = np.min(source.xmin)
                    zmax = np.max(source.xmax)
                ax = fig.add_subplot(grid[p%4, int(np.floor(p/4))])
                ax.set_ylabel("Flux [" + str(source._spec.y.unit) + "]")
                #if (hasattr(source._spec, '_orig')):
                #    spec = dc(source._spec._orig)
                #else:
                #    spec = dc(source._spec)
                spec = dc(source._spec)
                line = dc(source._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
                ax.plot(spec.x, spec.y, c='black', lw=1.0)
                #ax.plot(spec.x, spec.dy, c='r', lw=1.0)
                #ax.plot(spec.x, -spec.dy, c='r', lw=1.0)
                if (chunk is not None):
                    if (hasattr(source, '_norm')):
                        norm = dc(source._norm)
                        norm.to_z([ion[p]])
                        #for c in range(1, len(chunk)):
                        #    ax.plot(norm.x[chunk[c]], norm.y[chunk[c]], c='y',
                        #            lw=1.0, linestyle=':')
                        ax.plot(norm.x[chunk_sum], norm.y[chunk_sum], c='y',
                                lw=1.0, linestyle=':')
                    if (hasattr(source, '_voigt')):
                        voigt = dc(source._voigt)
                        voigt.to_z([ion[p]])
                        #for c in range(1, len(chunk)):
                        #    ax.plot(voigt.x[chunk[c]], voigt.y[chunk[c]], c='g',
                        #            lw=1.0, linestyle=':')
                        ax.plot(voigt.x[chunk_sum], voigt.y[chunk_sum], c='g',
                                lw=1.0, linestyle=':')
                    if (hasattr(source, '_cont')):
                        cont = dc(source._cont)
                        cont.to_z([ion[p]])
                        #for c in range(1, len(chunk)):
                        #    ax.plot(cont.x[chunk[c]], cont.y[chunk[c]], c='y')
                        ax.plot(cont.x[chunk_sum], cont.y[chunk_sum], c='y')
                    if (hasattr(source, '_fit')):
                        fit = dc(source._fit)
                        fit.to_z([ion[p]])
                        #for c in range(1, len(chunk)):
                        #    ax.plot(fit.x[chunk[c]], fit.y[chunk[c]], c='g')
                        ax.plot(fit.x[chunk_sum], fit.y[chunk_sum], c='g')
                    #if (hasattr(source, '_resid_fit')):
                    #    resid_fit = dc(source._resid_fit)
                    #    resid_fit.to_z([ion[p]])
                    #    ax.plot(resid_fit.x[chunk_sum],
                    #            resid_fit.y[chunk_sum], c='b', lw=1.0)
                    #if (hasattr(source, '_rem')):
                    #    rem = dc(source._rem)
                    #    rem.to_z([ion[p]])
                    #    ax.plot(rem.x[chunk_sum], rem.y[chunk_sum], c='b',
                    #            lw=1.0)
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
            if (hasattr(source, '_fit')):
                fig.suptitle("Reduced chi-squared: %3.2f" % (source._redchi),
                             fontsize=10)
        elif mode == 'compare':
            if (chunk is not None):
                zmin = np.min(source.xmin[group[1]])
                zmax = np.max(source.xmax[group[1]])
            else:
                zmin = np.min(source.xmin)
                zmax = np.max(source.xmax)
            grid = gs(1,1)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            ax = fig.add_subplot(grid[:,:])
            ax.set_xlabel("Redshift")
            ax.set_ylabel("Flux [" + str(source._spec.y.unit) + "]")
            for p in range(n):
                spec = dc(source._spec)
                line = dc(source._line)
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
            zmin = np.min(source.xmin)
            zmax = np.max(source.xmax)
            grid = gs(1,1)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            ax = fig.add_subplot(grid[:,:])
            ax.set_xlabel("Redshift")
            ax.set_ylabel("Flux [" + str(source._spec.y.unit) + "]")
            spec = dc(source._spec)
            line = dc(source._line)
            spec.to_z([ion[0]])
            line.to_z([ion[0]])
            ax.set_xlim(zmin, zmax)
            ax.plot(spec.x, spec.y, c='black', lw=1.0)
            ax.plot(spec.x, spec.dy, c='r', lw=1.0)
            if (chunk is not None):
                if (hasattr(source, '_norm')):
                    norm = dc(source._norm)
                    norm.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(norm.x[chunk[c]], norm.y[chunk[c]], c='y',
                                lw=1.0, linestyle=':')
                if (hasattr(source, '_voigt')):
                    voigt = dc(source._voigt)
                    voigt.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(voigt.x[chunk[c]], voigt.y[chunk[c]], c='g',
                                lw=1.0, linestyle=':')
                if (hasattr(source, '_cont')):
                    cont = dc(source._cont)
                    cont.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(cont.x[chunk[c]], cont.y[chunk[c]], c='y')
                if (hasattr(source, '_fit')):
                    fit = dc(source._fit)
                    fit.to_z([ion[0]])
                    for c in range(1, len(chunk)):
                        ax.plot(fit.x[chunk[c]], fit.y[chunk[c]], c='g')
                if (hasattr(source, '_rem')):
                    rem = dc(source._rem)
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
    """

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
        source = self #._source
        #print(self.t)
        #print(self._source.t)
        #print(self._flat.t)
        
        sumlen = len(z) + len(N) + len(b) + len(btur)
        if ((z != []) and (sumlen % len(z) != 0)):
            raise Exception("Parameter arrays must have the same length.")

        syst = dc(source)
        syst._t = source.t_all
        model = Model(self._spec, syst=syst, group=group, chunk=chunk)
        if (z == []):
            #z = self.x[group[1]]
            z = self._z[group[1]]
            if (hasattr(self, '_norm')):
                N = model.N_guess(self._norm, ion=source._flat.ion)
            else:
                N = model.N_guess(self._unabs, ion=source._flat.ion)
            if (hasattr(self, '_unabs')):
                cont = self._unabs
            elif (hasattr(self, '_norm')):
                cont = self._norm
            else:
                raise Exception("Continuum not found.")
            N = model.N_guess(cont, ion=source._flat.ion)
            b = np.full(len(source._t_all['X'][group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(source._t_all['X'][group[1]]), voigt_def['btur']) \
                   * u.km / u.s

        else:
            for val in N:
                if (val == voigt_def['N']):
                    N = model.N_guess(self._norm, ion=source._flat.ion)
            
        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)

        ion = np.unique(source._flat.ion)
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
