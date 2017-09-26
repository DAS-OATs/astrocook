from . import Line, Model
from .utils import convolve, dict_wave, voigt_def
from astropy import units as u
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
                 yunit=None,
                 meta=None,
                 dtype=float):
        """ Constructor for the Syst class """ 
        
        # Exceptions
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
            if (hasattr(line, '_cont')):
                self._precont = dc(line._precont)
                self._cont = dc(line._cont)            
            self._minima = dc(line._minima)
            self._maxima = dc(line._maxima)            

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

    def corr_resid(self, group, cont_corr):
        """ Add a new line at the minimum residual """

        neb = 'neb'
        
        zmin_ion = np.min(self.xmin[group[1]])
        zmax_ion = np.max(self.xmax[group[1]])
        where = self._chunk_sum
        resid_norm = self._resid_fit.y.value * 0.0
        resid_norm[where] = self._resid_fit.y[where]/self._resid_fit.dy[where]
        x = self._resid_fit.x[np.argmin(resid_norm)]
        ion_arr = np.unique(self._flat.ion)
        n = len(ion_arr)
        z_arr = np.empty(n)        
        xmin_ion = float('inf') * u.nm
        xmax_ion = 0 * u.nm
        for p in range(n):
            z_arr[p] = x / dict_wave[ion_arr[p]] - 1.0
            xmin_ion = min(xmin_ion, (1 + zmin_ion) * dict_wave[ion_arr[p]])
            xmax_ion = max(xmax_ion, (1 + zmax_ion) * dict_wave[ion_arr[p]])
        #print(xmin_ion, xmax_ion)
        where = abs(z_arr-np.mean([zmin_ion, zmax_ion])) \
                == abs(z_arr-np.mean([zmin_ion, zmax_ion])).min()
        z_ion = z_arr[where][0]
        z_neb = x / dict_wave[neb] - 1.0
        zmin_neb = xmin_ion / dict_wave[neb] - 1.0
        zmax_neb = xmax_ion / dict_wave[neb] - 1.0
        print(xmin_ion, xmax_ion)
        print(zmin_ion, zmax_ion)
        
        where = abs(z_ion-self.x.value) == abs(z_ion-self.x.value).min()
        ion = self.ion[where][0]
        y_ion = np.empty(len(ion))
        dy_ion = np.empty(len(ion))        
        for i in range(len(ion)):
            spec = dc(self._spec)
            spec.to_z([ion[i]])
            y_ion[i] = np.interp(z_ion, spec.x, spec.y)
            dy_ion[i] = np.interp(z_ion, spec.x, spec.dy)
        spec = dc(self._spec)
        spec.to_z([neb])
        y_neb = np.interp(z_neb, spec.x, spec.y)
        dy_neb = np.interp(z_neb, spec.x, spec.dy)
        xmin_neb = x

        #self._t.add_row([z_ion, y_ion, zmin_ion, zmax_ion, dy_ion, ion])
        self.flatten_z()
        #print(z_neb, y_neb, zmin_neb, zmax_neb, dy_neb, neb)
        self._flat.t.add_row([z_neb, y_neb, zmin_neb, zmax_neb, dy_neb, neb])
        #print(self._flat.t)
        self.deflatten_z()
        print(self._t)
        #sys.exit()
        
        self.t.sort('X')  # This gives an annoying warning

        return cont_corr
        
    def chunk(self, x=None, line=None):  # Chunk must be shifted to the system z
        if ((x is None) and (line is None)):
            raise Exception("Either x or line must be provided.")
        if (x is not None):
            if (line is not None):
                warnings.warn("x will be used; line will be disregarded.")
            line = np.where(abs(self.x-x.value) \
                            == abs(self.x-x.value).min())[0][0]
        if ((x is None) and (line >= len(self._t))):
            raise Exception("Line number is too large.")
        try:  # When ION has different sizes in different rows
            ion = np.unique(np.asarray(np.sum(self.ion)))
        except:
            ion = np.unique(self.ion)
        n = len(ion)
        iter = range(len(self._t))
        ret = (line,)
        #print(self.t)
        for p in range(n):
            sel = self._spec.t['X'] < 0.0
            spec = dc(self._spec)
            spec.to_z([ion[p]])
            #print(self.t[self.group(line=line)[1]])
            for row in self.t[self.group(line=line)[1]]:
                #print(row['XMIN'], row['XMAX'])
                sel = np.logical_or(sel, np.logical_and(
                    spec.t['X'] >= row['XMIN'],
                    spec.t['X'] <= row['XMAX']))
            if (np.sum(sel) % 2 == 0):
                sel[np.argmax(sel)] = 0

            #print(ion[p], np.sum(sel))
                
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

        z = 0
        end_row = False
        for l in range(len(self._flat.t)):
            if (np.isclose(self._flat.x[l], z) == False):
                if (end_row == True):
                    z_deflat.append(z)
                    y_deflat.append(y)
                    zmin_deflat.append(zmin)
                    zmax_deflat.append(zmax)
                    dy_deflat.append(dy)
                    ion_deflat.append(ion)
                    end_row = False
                y = (self._flat.y[l].value,)
                zmin = self._flat.xmin[l]
                zmax = self._flat.xmax[l]
                dy = (self._flat.dy[l].value,)                
                ion = (self._flat.ion[l],)
            else:
                z = self._flat.x[l]
                y = y + (self._flat.y[l].value,)
                dy = dy + (self._flat.dy[l].value,)
                ion = ion + (self._flat.ion[l],)
            z = self._flat.x[l]
            end_row = True
            #print(self._flat.x[l])
            #print(z_deflat)
        z_deflat.append(z)
        y_deflat.append(y)
        zmin_deflat.append(zmin)
        zmax_deflat.append(zmax)
        dy_deflat.append(dy)
        ion_deflat.append(ion)
        
        syst = Syst(self.line, self.spec, x=z_deflat, y=y_deflat,
                    xmin=zmin_deflat, xmax=zmax_deflat, dy=dy_deflat,
                    ion=ion_deflat, yunit=yunit)
        self.__dict__.update(syst.__dict__)
    
    def fit(self, group, chunk, unabs_guess, voigt_guess, psf, maxfev=1000):

        """
        for c in range(1, len(chunk)):
            #print(np.sum(chunk[c]))
            if (c == 1):
                #model = unabs_guess[2*(c-1)] * voigt_guess[2*(c-1)]
                model = model * voigt_guess[2*(c-1)]
                conv_model = lmc(model, psf[2*(c-1)], convolve)
                #param = unabs_guess[2*c-1]
                param.update(voigt_guess[2*c-1])
                #param = voigt_guess[2*c-1]
                param.update(psf[2*c-1])
                chunk_sum = dc(chunk[c])
            else:
                #model = unabs_guess[2*(c-1)] * voigt_guess[2*(c-1)]
                model = voigt_guess[2*(c-1)]
                #conv_model += lmc(model, psf[2*(c-1)], convolve)
                conv_model *= lmc(model, psf[2*(c-1)], convolve)
                #param.update(unabs_guess[2*c-1])
                param.update(voigt_guess[2*c-1])
                param.update(psf[2*c-1])
                chunk_sum += chunk[c]
        """
        model = unabs_guess[0] * voigt_guess[0]
        param = unabs_guess[1]
        param.update(voigt_guess[1])
        for c in range(1, len(chunk)):
        #    model *= voigt_guess[2*(c-1)]
        #    param.update(voigt_guess[2*c-1]) 
        #    #param.pretty_print()
            if (c == 1):
                chunk_sum = dc(chunk[c])
            else:
                chunk_sum += chunk[c]
        #conv_model = lmc(model, psf[0], convolve)
        #param.update(psf[1])
        conv_model = model
         
        #expr_dict = voigt_guess[-1]
        #for k in expr_dict:
        #    param[k].set(expr=expr_dict[k])

        if (hasattr(self, '_cont') == False):
            self._cont = dc(self._spec)            
        if (hasattr(self, '_fit') == False):
            self._fit = dc(self._spec)
        if (len(self._spec.x[chunk_sum]) < len(param)):
            warnings.warn("Too few data points; skipping.")
            fit = None
        else:
            if (hasattr(self, '_norm')):
                #print("cont")
                y = self._spec.y[chunk_sum] / self._cont.y[chunk_sum]
                dy = self._spec.dy[chunk_sum].value \
                     / self._cont.y[chunk_sum].value
                param.pretty_print()
                fit = conv_model.fit(y.value, param,
                                     x=self._spec.x[chunk_sum].value,
                                     fit_kws={'maxfev': maxfev},
                                     weights=1/dy)
                print(fit.fit_report())
                #print('tot', np.sum(chunk_sum))
                comp = fit.eval_components(x=self._spec.x[chunk_sum].value)
                self._fit.y[chunk_sum] = fit.best_fit \
                                         * self._cont.y[chunk_sum] \
                                         * self._fit.y[chunk_sum].unit
                #self._cont.y[chunk_sum] = comp['cont1_'] \
                #                          * self._cont.y[chunk_sum].value
                cont_temp = dc(self._cont)
                for c in range(1, len(chunk)):
                    #print(comp['cont' + str(c) + '_'])
                    """
                    if (c == 1):
                        #cont_temp.y[chunk_sum] = comp['cont' + str(c) + '_']
                        cont_temp.y[chunk_sum] = comp['cont_'] \
                                                 * self._cont.y[chunk_sum].value
                    else:
                        #cont_temp.y[chunk_sum] += comp['cont' + str(c) + '_'] 
                        cont_temp.y[chunk_sum] += comp['cont_'] \
                                                 * self._cont.y[chunk_sum].value
                    """
                    #print(c, np.mean(cont_temp.y[chunk_sum]))
                cont_temp.y[chunk_sum] = comp['cont_'] \
                                         * self._cont.y[chunk_sum].value
                self._cont = cont_temp
                
            else:
                fit = conv_model.fit(self._spec.y[chunk_sum].value, param,
                                     x=self._spec.x[chunk_sum].value,
                                     fit_kws={'maxfev': maxfev},
                                     weights=1/self._spec.dy[chunk_sum].value)

                cont = fit.eval_components(x=self._spec.x[chunk_sum].value)
                """
                for c in range(1, len(chunk)):
                    if (c == 1):
                        self._cont.y[chunk_sum] = cont['cont' + str(c) + '_'] \
                                                  * self._cont.y[chunk_sum].unit
                    else:
                        self._cont.y[chunk_sum] += cont['cont' + str(c) + '_'] \
                                                * self._cont.y[chunk_sum].unit
                """
                self._cont.y[chunk_sum] = cont['cont' + str(c) + '_'] \
                                          * self._cont.y[chunk_sum].unit
                self._fit.y[chunk_sum] = fit.best_fit \
                                         * self._fit.y[chunk_sum].unit
            self._redchi = fit.redchi
            self._aic = fit.aic
            self._chunk_sum = chunk_sum
        
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

        # This part is needed to add interlopers to the group
        xmin = float('inf') 
        xmax = 0
        for r in self._t[sel]:
            for i in range(len(r['Y'])):
                xmin = min(xmin, (1 + r['XMIN']) * dict_wave[r['ION'][i]])
                xmax = max(xmax, (1 + r['XMAX']) * dict_wave[r['ION'][i]])
        l = 0
        for r in self._t:
            for i in range(len(r['Y'])):
                x = (1 + r['X']) * dict_wave[r['ION'][i]]
                if ((x > xmin) and (x < xmax)):
                    sel[l] = True
            l += 1
        return (line, sel)
        
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
            * self._cont.y[chunk_sum] * self._norm.y[chunk_sum].unit 

        #print(np.mean(self._norm.y[chunk[1]]))
        return norm 
    def plot(self, group=None, chunk=None, figsize=(10,4), split=False,
             block=True, **kwargs):
        ion = np.unique(self._flat.ion)
        #ion = ion[ion != 'neb'] 
        n = len(ion)
        z = self.x
        if (chunk is not None):
            zmin = np.min(self.xmin[group[1]])
            zmax = np.max(self.xmax[group[1]])
        else:
            zmin = np.min(self.xmin)
            zmax = np.max(self.xmax)
        if split == True:
            #figsize = (7,6)
            row = min(n,4)
            col = int(np.ceil(n/4))
            figsize = (col*6, n*2.5)
            fig = plt.figure(figsize=figsize)
            fig.canvas.set_window_title("System")
            grid = gs(row,col)
            for p in range(n):
                t = Table(self._t[group[1]])
                #print(ion[p], t['ION'])
                for l in range(len(t)):
                    if (ion[p] in t['ION'][l]):
                        zmin = t['XMIN'][l]
                        zmax = t['XMAX'][l]
                ax = fig.add_subplot(grid[p%4, int(np.floor(p/4))])
                ax.set_ylabel("Flux [" + str(self._spec.y.unit) + "]")
                spec = dc(self._spec)
                line = dc(self._line)
                spec.to_z([ion[p]])
                line.to_z([ion[p]])
                ax.plot(spec.x, spec.y, c='black', lw=1.0)
                ax.plot(spec.x, spec.dy, c='r', lw=1.0)
                #if (hasattr(self, '_unabs')):
                #    unabs = dc(self._unabs)
                #    unabs.to_z([ion[p]])
                #    for c in range(1, len(chunk)):
                #        ax.plot(unabs.x[chunk[c]], unabs.y[chunk[c]], c='y',
                #                lw=1.0, linestyle=':')
                if (hasattr(self, '_norm')):
                    norm = dc(self._norm)
                    norm.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        ax.plot(norm.x[chunk[c]], norm.y[chunk[c]], c='y',
                                lw=1.0, linestyle=':')
                if (hasattr(self, '_voigt')):
                    voigt = dc(self._voigt)
                    voigt.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        ax.plot(voigt.x[chunk[c]], voigt.y[chunk[c]], c='g',
                                lw=1.0, linestyle=':')
                if (hasattr(self, '_cont')):
                    cont = dc(self._cont)
                    cont.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        #print(c, np.mean(cont.y[chunk[c]]))
                        ax.plot(cont.x[chunk[c]], cont.y[chunk[c]], c='y')
                if (hasattr(self, '_fit')):
                    fit = dc(self._fit)
                    fit.to_z([ion[p]])
                    for c in range(1, len(chunk)):
                        ax.plot(fit.x[chunk[c]], fit.y[chunk[c]], c='g')
                    ax.plot(fit.x, fit.y, c='g')
                ax.scatter(line.x, line.y, c='b')
                for comp in z:
                    ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='black')
                if ((p+1) % row != 0):
                    pass
                    #ax.set_xticks([], [])
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
                ax.plot(spec.x, spec.y, lw=1.0)
                ax.scatter(line.x, line.y, marker='+')
            text = ', '.join(str(p) for p in ion)
            for comp in z:
                ax.axvline(x=comp, ymin=0.65, ymax=0.85, color='black', lw=3.0,
                           linestyle='--')

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
            #print(N)
            b = np.full(len(self.x[group[1]]), voigt_def['b']) * u.km / u.s
            btur = np.full(len(self.x[group[1]]), voigt_def['btur']) \
                   * u.km / u.s

        if (hasattr(self, '_voigt') == False):
            self._voigt = dc(self._spec)

        ion = np.unique(self._flat.ion)
        voigt = model.voigt(z, N, b, btur, ion)
        #print(voigt)
        for c in range(1, len(chunk)):
            #print(len(voigt))
            #print(c, voigt[2*(c-1)])
            #voigt[2*c-1].pretty_print()
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
        return voigt
