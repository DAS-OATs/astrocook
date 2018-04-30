from . import Cont, Line, Model, Spec1D
from .utils import *
from .model import voigt_params
from astropy import units as u
from astropy.io import fits as fits
from astropy.table import Column, Table, join, unique, vstack
from copy import deepcopy as dc
from lmfit import CompositeModel as lmc
from matplotlib.gridspec import GridSpec as gs
import matplotlib.pyplot as plt
import numpy as np
import warnings
import sys
import time

class System(Spec1D, Line, Cont):

    def __init__(self, spec=None, line=None, cont=None,
                 series=None, ion=None,
                 z=None, N=None, b=None, btur=None,  
                 dz=None, dN=None, db=None, dbtur=None,
                 vary=[True,True,True,False], expr=[None,None,None,None],
                 x=None, y=None, xmin=None, xmax=None, dy=None,  
                 Nunit=Nunit_def, bunit=bunit_def,
                 xunit=xunit_def, yunit=yunit_def, 
                 meta=None, dtype=float):  
        """ Constructor for the System class """ 

        
        # Spectrum
        if (spec != None):
            self._spec = dc(spec)
        else:
            self._spec = dc(line._spec)

        # "is not" works also with arrays
        if (z is not None):            
            if (btur is None):
                btur = 0.0 * bunit_def
            self._t = self.create_t(series, z, N, b, btur,
                                    dz, dN, db, dbtur, vary, expr,
                                    Nunit, bunit, dtype)
            print self._t
            self.create_line(xmin=xmin, xmax=xmax)
            
        # Line list
        else:
            if (line != None):
                self._line = dc(line)
            elif (x is not None):
                if (x.unit == None):
                    x = x * xunit_def
                if (y.unit == None):
                    y = y * yunit_def
                self._line = Line(spec=spec, x=x, y=y, xmin=xmin, xmax=xmax,
                                  dy=dy)
            if (ion != None):
                self.create_z(ion)


        # Continuum
        if (cont != None):
            self._cont = dc(cont)
        else:
            self._cont = dc(line._cont)

        """
        # Line list -- Deprecated!
        self._line = None
        if (line != None):
            self._line = dc(line)
            # Continuum 
            if (hasattr(line, '_cont')):
                self._precont = dc(line._precont)
                self._cont = dc(line._cont)
            if (hasattr(line, '_minima')):
                self._minima = dc(line._minima)
            if (hasattr(line, '_maxima')):
                self._maxima = dc(line._maxima)
        #else:
        #    self._line = Line(spec=spec, x=x, y=y, xmin=xmin, xmax=xmax, dy=dy)
            
        # Spectrum -- Deprecated!
        self._spec = None
        if (spec is not None):
            #if (hasattr(spec, '_orig')):
            #    self._spec = dc(spec._orig)
            #else:
            #    self._spec = dc(spec)
            self._spec = dc(spec)
            
        # Ion list
        if (line is not None):
            len_ion = len(line._t)
        if (x != []):
            print("here")
            len_ion = len(x)
        #print(x, x.value, x is not [])

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
        if isinstance(value, System):
            self._linez = value
        else:
            raise Exception("Redshift list has a wrong format.")

    """
    @property
    def t(self):
        if self._use_good:
            return self._t[self._igood]
        else:
            return self._t        
    """
    
# Methods

    def add(self, series, z, **kwargs):
        """ Add a system to a list """

        add = self.create_t(series=series, z=z, **kwargs)
        self._t = vstack([self._t, add])
        #self.create_line()
        #self.group(z)
        #self.chunk(z)
        #self.model(z)
        
    def chunk(self, z, dx=0.0, **kwargs):
        """ Extract the chunks of spectrum needed for fitting """

        if (hasattr(self, '_group') == False or True):
            self.group(z, dx, **kwargs)

        if (hasattr(self, '_chunk') == False or True):
            spec = dc(self._spec.t)
            spec.add_column(Column(self._cont.t['Y'], name='CONT'))
        else:
            spec = dc(self._chunk)

        x = spec['X']
        where = np.array([], dtype=int)
        for l in self._group:
            xmin = l['XMIN']-dx
            xmax = l['XMAX']+dx
            #cond_temp = np.logical_and(x>xmin, x<xmax)
            cond_temp = np.logical_and(x>=xmin, x<=xmax)
            where_temp = np.where(cond_temp)[0]
            if (len(where_temp) % 2 == 0):
                where_temp = where_temp[:-1]
            where = np.append(where, where_temp)
        where = np.unique(where)
        self._chunk = spec[where]
        self._chunk.sort('X')

        # Table with indexes of the chunk limits
        imin = np.append(0, np.where(np.ediff1d(where)>1)[0]+1)
        imax = np.append(np.where(np.ediff1d(where)>1)[0]+1, len(self._chunk))
        imean = (imax+imin)//2
        self._chunk_lim = Table()
        self._chunk_lim['MIN'] = Column(imin, dtype=int)
        self._chunk_lim['MAX'] = Column(imax, dtype=int)
        self._chunk_lim['MEAN'] = Column(imean, dtype=int)
        
    def create_line(self, xmin=None, xmax=None, sigma=0.07):
        """ Create a list of lines from a list of systems """

        x_temp = np.array([])
        xmin_temp = np.array([])
        xmax_temp = np.array([])
        z_temp = np.array([])

        # Find adjacent maxima to determine xmin and xmax - Improve!
        conv = self._spec.convolve(gauss_sigma=sigma, convert=False)
        mins, maxs, exts = conv.find_extrema()
        self._conv = conv
        for s in self._t:
            z = s['Z']
            for i in s['SERIES']:
                wave = dict_wave[i]
                wave_z = wave.to(xunit_def) * (1+z)
                x_temp = np.append(x_temp, wave_z)
                z_temp = np.append(z_temp, z)
                if (xmin == None or xmax == None):
                    pos_min, wave_min = find_nearest(np.array(maxs['X']),
                                                     wave_z.value)
                    if wave_min > wave_z.value:
                        pos_min = pos_min-1
                    xmin_temp = np.append(xmin_temp, maxs['X'][pos_min])
                    xmax_temp = np.append(xmax_temp, maxs['X'][pos_min+1])
                
        #x = np.sort(x_temp)
        x = x_temp
        if (xmin == None or xmax == None):
            xmin = xmin_temp
            xmax = xmax_temp
        y = np.interp(x, self._spec.x, self._spec.y.to(yunit_def))
        dy = np.interp(x, self._spec.x, self._spec.dy.to(yunit_def))
        
        self._line = Line(spec=self._spec, x=x, y=y,
                          xmin=xmin, xmax=xmax,
                          dy=dy, xunit=xunit_def, yunit=yunit_def)
        
        # Table to map rows of self._line.t into rows of self._t
        self._map = Table()
        self._map['X'] = Column(x_temp.value, dtype=float, unit=xunit_def)
        self._map['Z'] = Column(z_temp, dtype=float, unit=u.nm/u.nm)

    def create_t(self, series=[['unknown']],
                 z=z_def, N=N_def, b=b_def, btur=btur_def, 
                 dz=None, dN=None, db=None, dbtur=None,
                 vary=[True, True, True, False], expr=[None,None,None,None],
                 Nunit=Nunit_def, bunit=bunit_def, dtype=float):
        """ Create a list of systems """

        z = np.array(z, ndmin=1)
        N = np.array(N, ndmin=1)
        b = np.array(b, ndmin=1)
        btur = np.array(btur, ndmin=1)
        dz = np.array(dz, ndmin=1)
        dN = np.array(dN, ndmin=1)
        db = np.array(db, ndmin=1)
        dbtur = np.array(dbtur, ndmin=1)
        vary = np.array(vary, ndmin=2)
        expr = np.array(expr, ndmin=2)
        zunit = u.nm/u.nm
        t = Table()
        t['Z'] = Column(z, dtype=dtype, unit=zunit)
        t['N'] = Column(N, dtype=dtype, unit=Nunit)
        t['B'] = Column(b, dtype=dtype, unit=bunit)
        t['BTUR'] = Column(btur, dtype=dtype, unit=bunit)
        t['DZ'] = Column(dz, dtype=dtype, unit=zunit)
        t['DN'] = Column(dN, dtype=dtype, unit=Nunit)
        t['DB'] = Column(db, dtype=dtype, unit=bunit)
        t['DBTUR'] = Column(dbtur, dtype=dtype, unit=bunit)
        t['VARY'] = Column(vary)
        t['EXPR'] = Column(expr)
        
        # Needed to have the series column without defined shape
        t.add_column(Column(dtype=object, length=len(t),
                            shape=1, name='SERIES'), index=0)
        for (i, s) in enumerate(t):
            s['SERIES'] = series[0]
            #s['SERIES'] = series[i]

        return t
        
    def create_z(self, ion):
        """ Create a list of redshifts from a list of lines """

        ion_in = np.array(ion, ndmin=1)
        if (len(ion_in.shape) > 1):
            raise Exception("Ion must be a scalar or a 1-d array.") 

        x = np.array([])
        ion = np.array([])
        z = np.array([])
        for i in ion_in:
            x_temp = self._line.x
            wave = dict_wave[i]
            x = np.append(x, x_temp.to(xunit_def).value)
            ion = np.append(ion, np.full(len(self._line.t), i))
            z = np.append(z, (x_temp/wave).value - 1)
        
        # Redshift table
        self._z = Table()
        self._z['X'] = Column(x, dtype=float, unit=xunit_def)
        self._z['ION'] = Column(ion) 
        self._z['Z'] = Column(z, dtype=float, unit=u.nm/u.nm)

    def extract(self, row):
        sel = self._t[row]
        syst_sel = System(
            spec=self._spec, cont=self._cont,
            series=[sel['SERIES']],
            z=sel['Z'], N=sel['N'], b=sel['B'], btur=sel['BTUR'],
            dz=sel['DZ'], dN=sel['DN'], db=sel['DB'], dbtur=sel['DBTUR'],
            vary=sel['VARY'], expr=sel['EXPR'])
        return syst_sel
            
    def find(self, tag, ztol=1e-4):
        """ Find systems by matching redshifts """
        
        ion = dict_series[tag]
        self.create_z(ion)
        self._z.sort('Z')
        z_arr = self._z['Z']
        ion_arr = self._z['ION']
        match = np.isclose(z_arr[1:], z_arr[:-1], atol=ztol) 
        dec = np.core.defchararray.not_equal(ion_arr[1:], ion_arr[:-1])
        match = np.logical_and(match, dec)
        z_sel = np.mean([z_arr[1:][match], z_arr[:-1][match]], axis=0)
        self._t = self.create_t(series=[ion], z=z_sel)
        

    def fit(self, z, save=True, **kwargs):
        """ Fit the model on a system """

        if (hasattr(self, '_fun') == False or True):
            self.model(z, **kwargs)

        x = self._chunk['X']
        y = self._chunk['Y']
        dy = self._chunk['DY']
        cont = self._chunk['CONT']

        fun = self._fun
        par = self._par
        #par.pretty_print()
        fit = fun.fit(y/cont, par, x=x, weights=cont/dy)
        par = fit.params
        #par.pretty_print()
        model = fit.eval(par, x=x) * cont
        #self._model = fit.eval(par, x=self._spec.t['X']) * self._spec.t['CONT']
        #model = fit.best_fit * cont
        
        self._fun = fun
        self._par = par
        self._fit = fit
        for l in self._group:
            pref = l['PREF']
            l['Z'] = par[pref+'_z'].value
            l['N'] = par[pref+'_N'].value
            l['B'] = par[pref+'_b'].value
            l['BTUR'] = par[pref+'_btur'].value
            l['DZ'] = par[pref+'_z'].stderr
            l['DN'] = par[pref+'_N'].stderr
            l['DB'] = par[pref+'_b'].stderr
            l['DBTUR'] = par[pref+'_btur'].stderr
        self._chunk['MODEL'] = model
        self._t = unique(self._group[self._t.colnames], keys='Z')

    def group(self, z, dx, **kwargs):
        """ Extract the lines needed for fitting (both from systems and not)
        and setup the fitting parameters """

        if (hasattr(self, '_map') == False or True):
            self.create_line(dx, **kwargs)

        self._line.t.sort('X')
        

        # Join systems and lines
        join_t = join(join(self._t, self._map), self._line.t)
        
        # Select lines at the system redshift
        join_z = join_t['Z']
        cond_z = z==join_z
        
        # Select lines close to the previously selected ones
        join_xmin = join_t['XMIN']
        join_xmax = join_t['XMAX']
        cond_x = np.full(len(join_t), False)
        for j in join_t[cond_z]:
            xmin = j['XMIN']
            xmax = j['XMAX']
            #cond_x += np.logical_and(join_xmax>xmin, join_xmin<xmax)
            cond_x += np.logical_and(join_xmax>=xmin, join_xmin<=xmax)
        
        group = join_t[np.where(cond_x)[0]]
        group.sort(['Z', 'X'])

        # Define the ion and prefix columns
        ion = np.array([])
        for g in group:
            xs = (1+g['Z'])*np.array([dict_wave[i].value for i in g['SERIES']])
            ion = np.append(ion, g['SERIES'][np.abs(xs - g['X']).argmin()])
        pref = np.array(['voigt_'+str(i) for i in range(len(group))])
        group.add_column(Column(ion, name='ION'), index=1)
        group.add_column(Column(pref, name='PREF'), index=2)
        zlist = np.array(group['Z'])
        
        # Find rows with duplicate redshift values
        diff1d = np.append(zlist[0], np.ediff1d(zlist))  
        where = np.where(diff1d == 0)[0]

        # Associate each duplicate row to its companion, to link parameters 
        for (l, w) in enumerate(where):

            # Change the expression used in fit
            p = group['PREF'][w-1]
            group['VARY'][w] = [False, False, False, False]
            group['EXPR'][w] = [p+'_z', p+'_N', p+'_b', p+'_btur']

        self._group = group

    def model(self, z, norm=True, prof=True, psf=True, **kwargs):
        """ Create a model to fit, including normalization, profile and PSF """

        if (hasattr(self, '_chunk') == False or True):
            self.chunk(z, **kwargs)
            
        x = self._chunk['X']
        y = self._chunk['Y']
        dy = self._chunk['DY']
        cont = self._chunk['CONT']

        mods = Model(self._spec, syst=self)

        if (norm == True):
            mods.norm_new()#vary=[False])
            fun = mods._norm_fun
            par = mods._norm_par
        else:
            mods.norm_new(vary=[False])
            fun = mods._norm_fun
            par = mods._norm_par
            
        if (prof == True):
            # The only available profile shape is Voigt
            for (i, l) in enumerate(self._group):
                mods.voigt_new(ion=l['ION'],
                               z=l['Z'], N=l['N'], b=l['B'], btur=l['BTUR'],
                               vary=l['VARY'], expr=l['EXPR'], pref=l['PREF'])
                fun *= mods._prof_fun
                par.update(mods._prof_par)

        if (psf == True):
            for i, c in enumerate(self._chunk_lim):
                c_min = c['MIN']
                c_max = c['MAX']
                c_mean = c['MEAN']
                center = self._chunk['X'][c_mean]
                resol = self._chunk['RESOL'][c_mean]
                #print(center, resol)
                mods.psf_new2(c_min, c_max, center, resol, vary=False,
                              pref='psf_'+str(i))
                if i == 0:
                    psf_fun = mods._psf_fun
                    psf_par = mods._psf_par
                else:
                    psf_fun += mods._psf_fun
                    psf_par.update(mods._psf_par)

            #psf_y = np.concatenate(psf_fun.eval(psf_par, x=x))
            fun = lmc(fun, psf_fun, conv)
            par.update(psf_par)

        model = fun.eval(par, x=x) * cont
        
        try:
            self._chunk.add_column(Column(model, name='MODEL', dtype=float))
        except:
            self._chunk['MODEL'] = model
        self._fun = fun
        self._par = par
        
    def plot(self, ax=None, dz=0.01):
        """ Plot a system """

        if (hasattr(self, '_chunk') == False):
            self.chunk(z, **kwargs)

        t = self._t
        series = t['SERIES']
        ions = np.unique([[i for i in s] for s in series])
        waves = [dict_wave[i].value for i in ions]
        ions = ions[np.argsort(waves)]
        if ax == None:
            rown = 5
            n = len(ions)
            row = min(n,rown)
            col = int(np.ceil(n/rown))
            fig = plt.figure(figsize=(col*6, n*3.5))
            grid = gs(row,col)
            ax_arr = []
            axt_arr = []
            top = []
            bottom = []
            left = []
            for p in range(n):
                if p%rown == 0:
                    top.append(p)
                if p%rown==rown-1 or p==n-1:
                    bottom.append(p)
                if n<rown:
                    left.append(p)
                else:
                    left = [rown//2]
                ax_arr.append(fig.add_subplot(grid[p%rown,
                                                   int(np.floor(p/rown))]))
                axt_arr.append(ax_arr[p].twiny())
            try:
                fig.suptitle(r"$\chi_r^2$ = %3.1f" % self._fit.redchi)
            except:
                pass
                
        spec = dc(self._spec.t)
        group = dc(self._group)
        chunk = dc(self._chunk)
        x = spec['X']
        y = spec['Y']
        zmin = np.mean(group['Z'])-dz
        zmax = np.mean(group['Z'])+dz
        for c, i in enumerate(ions):
            x_z = x/dict_wave[i] - 1
            xmin = (1+zmin) * dict_wave[i].value
            xmax = (1+zmax) * dict_wave[i].value

            where_g = group['ION']==i
            z = group['Z'][where_g]

            where_c = np.where(np.logical_and(chunk['X']>=xmin,
                                              chunk['X']<=xmax))
            xc = chunk['X'][where_c]
            yc = chunk['Y'][where_c]
            contc = chunk['CONT'][where_c]
            modelc = chunk['MODEL'][where_c]
            residc = yc-modelc
            xc_z = xc / dict_wave[i] - 1
            
            ax_arr[c].set_xlim(zmin, zmax)
            axt_arr[c].set_xlim(xmin, xmax)
            maxf = 1.25
            ax_arr[c].set_ylim(-max(contc)*0.2, max(contc)*maxf)
            if c in bottom:
                ax_arr[c].set_xlabel(r"$z$")
            else:
                ax_arr[c].set_xticks([], [])
            if c in left:
                ax_arr[c].set_ylabel(r"Flux [%s]" % y.unit)
            axt_arr[c].tick_params(axis="x",direction="in", pad=-20)
            ax_arr[c].text(0.05, 0.5, i, transform=ax_arr[c].transAxes,
                           fontsize=13)
            ax_arr[c].plot(x_z, y, color='C0', linestyle='--')
            ax_arr[c].plot(x_z, self._conv.y, color='C2', linestyle=':')
            ax_arr[c].plot(xc_z, yc, color='C0')
            ax_arr[c].plot(xc_z, contc, color='C1')
            ax_arr[c].plot(xc_z, modelc, color='C2')
            ax_arr[c].plot(xc_z, residc, color='C3', linestyle=':')
            for g in group[where_g]:
                if c in top:
                    ax_arr[c].text(g['Z'], max(contc)*(maxf+0.18), 
                                   "%3.1f" % np.log10(g['N']), ha='center',
                                   fontsize=9) 
                    ax_arr[c].text(g['Z'], max(contc)*(maxf+0.1), 
                                   "%3.1f" % g['B'], ha='center', fontsize=9) 
                    ax_arr[c].text(g['Z'], max(contc)*(maxf+0.02), 
                                   "%3.1f" % g['BTUR'], ha='center', fontsize=9)
                ax_arr[c].axvline(x=g['Z'], color='C3', alpha=0.5)
        if ax == None:
            grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.9])
            grid.update(wspace=0.2, hspace=0.0)
            plt.show()

    def save(self, filename):
        #hdu = fits.BinTableHDU.from_columns(
        #    [self._t[i] for i in self._t.columns])
        #hdu.writeto(filename, overwrite=True)
        self._t.write(filename, format='fits', overwrite=True)
            
####
        
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
    """    
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
    """
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

    """
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
    """
        
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
                # Not working with Python 2.7
                #print("(%i) %3.2f;" \
                #      % (i, self._redchi), end=" ", flush=True)

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
        # Not working with Python 2.7
        #print("best chi-squared (%i) %3.2f, %3.2f; " \
        #      % (i_best, redchi_best, self._aic), end=" ", flush=True)

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
                # Not working with Python 2.7
                #print("(%i) %3.2f;" \
                #      % (i, self._redchi), end=" ", flush=True)
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
        # Not working with Python 2.7
        #print("best chi-squared (%i) %3.2f, %3.2f;" \
        #      % (i_best, redchi_best, self._aic), end=" ", flush=True)


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
            # Not working with Python 2.7
            #print("Redshift %i (%i/%i) (%3.4f)..." \
            #      % (l+1, l+1-list_range[0], len(list_range), x_arr[l].value),
            #      end=" ", flush=True)

            # Check if the group is new
            if (np.array_equal(self_temp.group(x=x_arr[l])[1], group_check)):
                print("same group, skipping.")
            else:
                group_check = self_temp.group(x=x_arr[l])[1]
                for i in iter_range:
                    
                    #self.fit_add(x=x_arr[l], i_max=i, mode=mode, **kwargs)
                    self.fit_auto(x=x_arr[l], i_max=i, mode=mode, **kwargs)
                    # Not working with Python 2.7
                    #print("time: %3.2f;" % (time.time()-start), end=" ",
                    #      flush=True)
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

        print(self._t)
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

    """
    def save(self, name):
        #hdu = fits.BinTableHDU.from_columns(
        #    [fits.Column(name='XMIN', format='E', array=self._spec.xmin),
        #     fits.Column(name='XMAX', format='E', array=self._spec.xmax),
        #     fits.Column(name='X', format='E', array=self._spec.x),
        #     fits.Column(name='Y', format='E', array=self._spec.y),
        #     fits.Column(name='Y_FIT', format='E', array=self._fit.y),
        #     fits.Column(name='Y_REM', format='E', array=self._rem.y),
        #     fits.Column(name='DY', format='E', array=self._spec.dy),
        #     fits.Column(name='GROUP', format='I', array=self._spec.group),
        #     fits.Column(name='RESOL', format='E', array=self._spec.resol)]) 
        #hdu.writeto(name + '_syst_spec.fits', overwrite=True)
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
    """
        
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
