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

    def __init__(self, acs=None, spec=None, line=None, cont=None,
                 series=None, ion=None,
                 z=None, N=N_def, b=b_def, btur=btur_def,  
                 dz=None, dN=None, db=None, dbtur=None,
                 vary=[True,True,True,False], expr=[None,None,None,None],
                 x=None, y=None, xmin=None, xmax=None, dy=None,  
                 Nunit=Nunit_def, bunit=bunit_def,
                 xunit=xunit_def, yunit=yunit_def, 
                 meta=None, dtype=float):  
        """ Constructor for the System class """ 

        if acs != None:
            self._acs(acs)

        
        # Spectrum
        if (spec != None):
            self._spec = dc(spec)
        #else:
        #    self._spec = dc(line._spec)


        # Continuum
        if (cont != None):
            self._cont = dc(cont)

            
        # "is not" works also with arrays
        if (z is not None):
            if (btur is None):
                btur = 0.0 * bunit_def
            self._t = self.create_t(series, z, N, b, btur,
                                    dz, dN, db, dbtur, vary, expr,
                                    Nunit, bunit, dtype)
            #if (self._line is None):
            #    self.create_line(xmin=xmin, xmax=xmax)
            
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


        self._use_good = False

    def _acs(self, acs):
        self.acs = acs
        self._spec = acs.spec
        self._line = acs.line
        self._cont = acs.cont
        self._model = acs.model
        #self._model._t['X'].mask = np.isnan(self._model._t['X'])
        try:
            self._map = acs.line._map  # When loading from a current session
        except:
            pass  # When loading from a saved session
        
        # Model and residuals
        """
        mask = dc(np.array(self._spec.t.mask))
        self._model = dc(self._spec.t)
        self._res = dc(self._spec.t)
        self._model.mask['X'] = True
        self._res.mask['X'] = True
        # This trick is needed to make the tables independent objects and
        # restore the original spectrum to unmasked state. Bug in Astropy?
        self._model['X'] = self._model['X']*-1*-1
        self._res['X'] = self._res['X']*-1*-1
        self._spec.t.mask = mask
        """
        
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

    def chunk(self): #, z, dx=0.0, **kwargs):
        """ @brief Extract a spectral chunk to fit a previously extracted group
        of lines
        """

        # Create a spectrum with attached continuum
        spec = dc(self.acs.spec)
        spec.t.add_column(Column(self.acs.cont.t['Y'], name='CONT'))

        # Regions around lines in the group are selected
        x = spec.t['X']
        where = np.array([], dtype=int)
        for l in self._group:
            xmin = l['XMIN']
            xmax = l['XMAX']
            cond_temp = np.logical_and(x>=xmin, x<=xmax)
            where_temp = np.where(cond_temp)[0]
            if (len(where_temp) % 2 == 0):
                where_temp = where_temp[:-1]
            where = np.append(where, where_temp)

        self._chunk_rows = np.unique(where)
        
        # Chunk is created as copy of the spectrum, completely masked apart
        # from the regions around lines
        self._chunk = spec.apply_mask(
            ['X', 'X'], [range(len(spec.t)), self._chunk_rows], [True, False]).t

        
    def extract_resid(self, s=None):
        """ @brief Extract residuals

        @param col Column used for masking regions where models is not defined
        """

        #out = dc(self.acs.model)
        out = self._model.apply_mask(
            ['X', 'X'], [range(len(self._model.t)), self._chunk_rows],
            [True, False])
        out.t['Y'] = out.t['YRESID']
        #out.t.remove_rows(out.t['X'].mask.nonzero()[0])
        out.t['Y'][out.t['X'].mask] = np.nan
        out.t.remove_columns(['YRESID', 'YADJ'])
        return out

        
    def fit(self, s=None, **kwargs):
        """ @brief Fit a model on a group of lines 
        
        @param s A row from a system table
        """

        if s == None:
            try:
                s = self._syst_sel
            except:
                s = self.acs._syst_sel
        
        self.model(s, **kwargs)

        #print self._group['EXPR', 'PREF']
                
        where = self._chunk['X'].mask == False
        x_c = self._chunk['X'][where]
        y_c = self._chunk['Y'][where]
        dy_c = self._chunk['DY'][where]
        cont_c = self._chunk['CONT'][where]

        fun = self._fun
        par = self._par
        #par.pretty_print()
        #print par['voigt_000_N'].value
        """
        print len(x_c)
        print len(y_c)
        print len(cont_c)
        print len(x_c[np.isnan(x_c)])
        print len(y_c[np.isnan(y_c)])
        print len(cont_c[np.isnan(cont_c)])
        """
        fit = fun.fit(np.array(y_c/cont_c), par, x=np.array(x_c))#,
                      #weights=np.array(cont_c/dy_c))
        #print fit.fit_report()
        par = fit.params
        #print par
        #print par['voigt_000_N'].value
        #print par['voigt_000_N'].stderr        
        #par.pretty_print()
        y = fit.eval(par, x=x_c) * cont_c
        yresid = y_c-y
        yadj = fit.eval_components()['adj_'] * cont_c
        #self._model.t['X'][:] = np.nan  # Otherwise they are converted into 1.0
                                        # when saving
        self._model.t['X'][where] = x_c
        self._model.t['Y'][where] = y
        self._model.t['DY'][where] = dy_c
        self._model.t['YRESID'][where] = yresid
        self._model.t['YADJ'][where] = yadj        
        self._fun = fun
        self._par = par
        self._fit = fit
        l_del = []
        for i, l in enumerate(self._group):
            #"""
            pref = l['PREF']
            if par[pref+'_z'].value > l['ZMIN'] and \
               par[pref+'_z'].value < l['ZMAX']:
                try:
                    l['Z'] = par[pref+'_z'].value
                    l['N'] = par[pref+'_N'].value
                    l['B'] = par[pref+'_b'].value
                    l['BTUR'] = par[pref+'_btur'].value
                    l['DZ'] = par[pref+'_z'].stderr
                    l['DN'] = par[pref+'_N'].stderr
                    l['DB'] = par[pref+'_b'].stderr
                    l['DBTUR'] = par[pref+'_btur'].stderr
                    l['X'] = (1.+l['Z'])*dict_wave[l['ION']].value
                except:
                    print "hey"
                    pass
            else:
                l_del.append(i)
        #self._group.remove_rows(l_del)
        #"""
            
        # Save fitted lines in the group
        try:
            cond = self._group['Z'] > 0
            new_t = unique(self._group[self._t.colnames][cond], keys='Z')
            #print new_t
            new_map = self._group[self._map.colnames][cond]
            new_line = self._group[self._line.t.colnames]
            self._t[self._group_t] = new_t#[new_t['Z']>0]
            self._map[self._group_map] = new_map
            self._line._t[self._group_line] = new_line
            self._line._t.sort('X')
        except:
            print "big hey"

        self.acs.model = self._model
        
    def group(self, s, **kwargs):
        """ @brief Crete a group of lines to be fitted together
        
        A group of line includes:
        1. Lines of the selected system;
        2. Lines of other systems with overlapping wavelength range;
        3. Lines of systems close in redshift with those at point 2.
        4. Other non-identified lines lines with overlapping wavelength range.
        A group of lines includes information about the guess fitting parameters
        of the lines and their constraints

        @param s A row from a system table
        """

        print np.array(s['Z'])
        
        self._line.t.sort('X')
        #print self._map

        # Join systems and lines
        join_t = join(join(self._t, self._map), self._line.t)
        
        # Select the system redshift 
        join_z = join_t['Z']
        cond_z = s['Z']==join_z
        #group = join_t[cond_z]
        
        # Select other systems close in redshift
        diff = 1
        while diff > 0:
            # First find the lines whose fitting ranges overlap with those of the
            # system lines..
            join_xmin = join_t['XMIN']
            join_xmax = join_t['XMAX']
            cond_x = np.full(len(join_t), False)
            for j in join_t[cond_z]:
                xmin = j['XMIN']
                xmax = j['XMAX']
                #cond_x += np.logical_and(join_xmax>=xmin, join_xmin<=xmax)
                cond_x += np.logical_and(join_xmax>xmin, join_xmin<xmax)

            diff = np.abs(np.sum(cond_x)-np.sum(cond_z))
            print diff
                
            # ...then select the whole systems those lines belong to...
            #join_zc = join_t['Z']
            cond_z = np.full(len(join_t), False)
            for j in join_t[np.where(cond_x)[0]]:
                z = j['Z']
                cond_z += z==join_z

        """
        # ...then select the lines whose fitting ranges overlap with those of
        # these systems...
        cond_xzc = np.full(len(join_t), False)
        for j in join_t[np.where(cond_zc)[0]]:
            xmin = j['XMIN']
            xmax = j['XMAX']
            cond_xzc += np.logical_and(join_xmax>xmin, join_xmin<xmax)

        # ...then select the whole systems those lines belong to
        cond_zxzc = np.full(len(join_t), False)
        for j in join_t[np.where(cond_xzc)[0]]:
            z = j['Z']
            print z
            cond_zxzc += z==join_z

        print np.sum(cond_x), np.sum(cond_zc), np.sum(cond_xzc), \
            np.sum(np.logical_or(cond_zc, cond_xzc)), \
            np.sum(cond_xzc)==np.sum(np.logical_or(cond_zc, cond_xzc))
        """
            
        group = join_t[cond_z]
        
        
        # Find wavelength duplicates (it may happen that a line may be
        # associated to two systems as the same ion) and remove them from group
        # and map
        group.sort('X')
        #print group
        #print np.where(np.ediff1d(group['X']) == 0.0)[0]
        group_where = np.where(np.ediff1d(group['X']) == 0.0)[0]
        map_where = np.where(np.logical_and(
            self._map['X'] == group['X'][group_where],
            self._map['Z'] == group['Z'][group_where]))
        group.remove_rows(group_where)
        self._map.remove_rows(map_where)

        # Sort by ascending redshift
        group.sort(['Z', 'X'])

        # Define the ion and prefix columns
        ion = np.array([])
        for ig, g in enumerate(group):
            series = dict_series[g['SERIES']]
            xs = (1+g['Z'])*np.array([dict_wave[i].value for i in series])
            ion = np.append(ion, series[np.abs(xs - g['X']).argmin()])
        pref = np.array(['voigt_%03d' % i for i in range(len(group))])
        group.add_column(Column(ion, name='ION'), index=1)
        group.add_column(Column(pref, name='PREF', format='20s'), index=2)
        zlist = np.array(group['Z'])

        # Define the mininimum- and maximum-redshift columns
        zmin = np.array([group['XMIN'][i]/dict_wave[group['ION'][i]].value-1
                         for i in range(len(group))])
        zmax = np.array([group['XMAX'][i]/dict_wave[group['ION'][i]].value-1
                         for i in range(len(group))])
        group.add_column(Column(zmin, name='ZMIN'), index=1)
        group.add_column(Column(zmax, name='ZMAX'), index=2)
    
        # Find rows with duplicate redshift values
        #print zlist
        diff1d = np.append(zlist[0], np.ediff1d(zlist))  
        where = np.where(diff1d == 0)[0]

        # Associate each duplicate row to its companion, to link parameters 
        for (l, w) in enumerate(where):

            # Change the expression used in fit
            p = group['PREF'][w-1]
            group['VARY'][w] = [False, False, False, False]
            group['EXPR'][w] = [p+'_z', p+'_N', p+'_b', p+'_btur']

        # Add unidentified lines close to lines in the system
        # Addition is iterative: also lines close to added lines are added,
        # until no new lines are found
        #iter_flag = True
        #sel_un = 0
        #while iter_flag:
        #"""
        for iter in []: #[0, 1]:
            #print "here"
            group_xmin = group['XMIN']
            group_xmax = group['XMAX']
            cond_un = np.full(len(self._line.t), False)
            for g in group:
                x = g['X']
                xmin = g['XMIN']
                xmax = g['XMAX']
                #cond_un += np.logical_and(self._line.t['XMAX']>=xmin,
                #                          self._line.t['XMIN']<=xmax)        
                cond_un += np.logical_and(self._line.t['XMAX']>xmin,
                                          self._line.t['XMIN']<xmax)
            i = len(group)
            for l in self._line.t[cond_un]:
                x = l['X']
                if np.all(x != group['X']):
                    lmap = self._map[self._map['X']==x]
                    xmin = l['XMIN']
                    xmax = l['XMAX']            
                    y = l['Y']
                    dy = l['DY']
                    ew = l['EW']
                    N = N_def
                    b = b_def
                    btur = btur_def
                    dz = None
                    dN = None
                    db = None
                    dbtur = None
                    vary = [True, True, True, True]
                    expr = [None, None, None, None]
                    if len(lmap) == 1:
                        z = lmap['Z']
                        series = self._t['SERIES'][self._t['Z']==z][0]
                        wave = np.array([dict_wave[s].value
                                         for s in dict_series[series]])
                        """
                        xs = (1+lmap['Z'])*wave
                        sel = np.abs(xs - lmap['X']).argmin()
                        xsel = wave[sel]
                        ion = dict_series[series][sel]
                        zmin = xmin/xsel-1
                        zmax = xmax/xsel-1
                        """
                        for c, s in enumerate(dict_series[series]):
                            pref = 'voigt_%03d' %i
                            wave_s = (1+lmap['Z']) * dict_wave[s].value
                            sel = np.abs(wave_s - self.line._t['X']).argmin()
                            x_s = self.line._t['X'][sel]
                            xmin_s = self.line._t['XMIN'][sel]
                            xmax_s = self.line._t['XMAX'][sel]
                            ion = s #dict_series[series][sel]
                            zmin = xmin_s/dict_wave[s].value-1
                            zmax = xmax_s/dict_wave[s].value-1
                            if c == 0:
                                p = pref
                            else:
                                vary = [False, False, False, False]
                                expr = [p+'_z', p+'_N', p+'_b', p+'_btur']
                            group.add_row([z, zmin, zmax, ion, pref, N, b, btur,
                                           dz, dN, db, dbtur, vary, expr,
                                           series, x_s, xmin_s, xmax_s, y, dy,
                                           ew])
                            i += 1
                    else:
                        pref = 'voigt_%03d' %i
                        series = 'unknown'
                        ion = 'unknown'
                        z = 0.0
                        zmin = xmin/x-1
                        zmax = xmax/x-1
                        group.add_row([z, zmin, zmax, ion, pref, N, b, btur, dz,
                                       dN, db, dbtur, vary, expr, series, x,
                                       xmin, xmax, y, dy, ew])
                        i += 1
                    """
                    group.add_row([z, zmin, zmax, ion, pref, N, b, btur, dz,
                                   dN, db, dbtur, vary, expr, series, x,
                                   xmin, xmax, y, dy, ew])
                    """
                    for s in group[-1:]:
                        self.N(s)
            #iter_flag = np.sum(cond_un) != sel_un
            #if np.sum(cond_un) == sel_un:
            #    break
            #sel_un = np.sum(cond_un)
        #"""
        self._group = group
        #print group
        
        self._t.sort('Z')
        self._map.sort('Z')
        self._group_t = np.in1d(self._t['Z'], group['Z'])
        self._group_map = np.in1d(self._map['Z'], group['Z'])
        self._group_line = np.in1d(self._line.t['X'], group['X'])
        #print self._t[self._group_t]
        #print self.line._t[self._group_line]

    def line_new(self, series='Ly_ab', mode='all'):
        """ @brief Use matching redshift from a list of lines to create a list
        of systems

        @param series Label of the series of transitions
        @param mode 'all': all line redshifts are defined as systems;
                    'match': only matching redshifts are defined as systems;
                    'complete': matching redshiftd are defined as systems, and
                    then remaining line redshifts are added to the system list
        """ 

        line = self.acs.line
        line_map = Table()
        if mode == 'all':
            z = np.array(line._z['Z'])
            line_map['X'] = line._z['X']
            line_map['Z'] = line._z['Z']
        if mode == 'match' or mode == 'complete':
            z = line._z_match
            line_map['X'] = np.append(line._z['X'][1:][line._w_match],
                                      line._z['X'][:-1][line._w_match])
            line_map['Z'] = np.append(line._z_match, line._z_match)
        out = System(acs=self.acs, series=series, z=z)
        line_map_add = Table()
        
        # If series is Ly_ab(...), discarded lines are added as Ly_a's
        if mode == 'complete':
            #print line._z['X']
            #print line_map['X']
            #print line._z
            #print dict_series[series][-1]
            disc = np.logical_and(~np.in1d(line._z['X'], line_map['X']),
                                  line._z['ION'] == dict_series[series][-1])
            #print line._z['X'][disc]
            #print line._z['X'][disc][0:460]
            z_add = line._z['Z'][disc]
            line_map_add['X'] = line._z['X'][disc]
            line_map_add['Z'] = line._z['Z'][disc]
            line_map = vstack([line_map, line_map_add])
            out_add = System(acs=self.acs, series=dict_series[series][-1],
                             z=z_add)
            out._t = vstack([out._t, out_add._t])
            out._t.sort('Z')
        line_map.sort('Z')
        out._map = line_map 
        #print out._t
        return out

    
    def model(self, s=None,
              adj='linear', adj_value=[1.0, 0.0], adj_vary=[False, False],
              adj_min=[None, None], adj_max=[None, None], adj_expr=[None, None],
              prof='voigt', psf='psf_gauss', psf_resol_value=-1.0,
              psf_resol_vary=False, psf_resol_min=None, psf_resol_max=None,
              psf_resol_expr=None, **kwargs):
        """ @brief Create a model for a group of lines, including continuum 
        adjustment, line profile and instrument PSF 

        @param s A row from a system table
        @param adj Type of continuum adjustment ('linear')
        @param adj_value Array of guess values for the parameters of adj
        @param adj_vary Array of constraint on variability for the parameters of
                        adj
        @param adj_min Array of minimum values for the parameters of adj
        @param adj_max Array of maximum values for the parameters of adj
        @param adj_expr Array of constraining expression for the parameters
        @param prof Type of line profile ('voigt')
        @param psf Type of instrument PSF ('psf_gauss')
        @param psf_resol_value Guess value for the PSF resolution
        @param psf_resol_vary Constraint on variability for the PSF resolution
        @param psf_resol_min Minimum value for the PSF resolution
        @param psf_resol_max Maximum value for the PSF resolution
        @param psf_resol_expr Constraining expression for the PSF resolution
        """

        if s == None:
            try:
                s = self._syst_sel
            except:
                s = self.acs._syst_sel
        
        # Determine group and chunk
        self.group(s)
        self.chunk()
        
        mod = Model()#self._chunk)

        # Continuum adjustment
        getattr(mod, adj)(adj_value, adj_vary, adj_min, adj_max, adj_expr)
        fun_adj = getattr(mod, '_'+adj+'_fun')
        par_adj = getattr(mod, '_'+adj+'_par')
        fun = fun_adj
        par = par_adj

        # Line profile
        for l in self._group:
            ion = l['ION']
            if ion == 'unknown':
                wave = np.array(l['X'])
            else:
                wave = 0.0
            prof_value = np.array([l['Z'], l['N'], l['B'], l['BTUR']])
            #print prof_value
            prof_vary = l['VARY']
            prof_expr = l['EXPR']
            prof_pref = l['PREF']
            getattr(mod, prof)(ion, wave, value=prof_value, vary=prof_vary,
                               expr=prof_expr, pref=prof_pref)        
            fun *= getattr(mod, '_'+prof+'_fun')
            par.update(getattr(mod, '_'+prof+'_par'))

        # Instrument PSF
        # The PSF cannot be convoluted on the masked spectrum; a spectrum with
        # only the non masked regions should be created, and the PSF should be
        # convoluted on each region at a time
        nz = self._chunk['X'].mask.nonzero()[0]
        cs = np.append(0, np.cumsum(nz[np.where(np.ediff1d(nz)>1)[0]+1]\
                                    -nz[np.where(np.ediff1d(nz)>1)[0]]-1))
        if len(cs) == 1:  # Case with a single chunk
            cs = np.array([0, len(self._chunk)-len(nz)])
        imin = cs[:-1]
        imax = cs[1:]
        imean = (imax+imin)//2
        rem = dc(self._chunk)
        rem.remove_rows(self._chunk['X'].mask.nonzero()[0])
        for i, (c_min, c_max, c_mean) in enumerate(zip(imin, imax, imean)):
            if psf_resol_value < 0.0:
                psf_resol_value = rem['RESOL'][c_mean]
            psf_value = [c_min, c_max, rem['X'][c_mean], psf_resol_value]
            psf_vary = [False, False, False, psf_resol_vary]
            psf_min = [None, None, None, psf_resol_min]
            psf_max = [None, None, None, psf_resol_max]
            psf_expr = [None, None, None, psf_resol_expr]
            psf_pref = 'psf_'+str(i) 
            getattr(mod, psf)(psf_value, psf_vary, psf_min, psf_max, psf_expr,
                               psf_pref)
            if i == 0:
                psf_fun = getattr(mod, '_'+psf+'_fun')
                psf_par = getattr(mod, '_'+psf+'_par')
            else:
                psf_fun += getattr(mod, '_'+psf+'_fun')
                psf_par.update(getattr(mod, '_'+psf+'_par'))

        fun = lmc(fun, psf_fun, conv)
        par.update(psf_par)

        # Update the model table
        where = self._chunk['X'].mask == False
        x_c = self._chunk['X'][where]
        y_c = self._chunk['Y'][where]
        dy_c = self._chunk['DY'][where]
        cont_c = self._chunk['CONT'][where]
        if len(self._model._t) == 0:
            x_t = self._chunk['X']
            xmin_t = self._chunk['XMIN']
            xmax_t = self._chunk['XMAX']
            y_t = np.empty(len(x_t))
            y_t[:] = np.nan
            dy_t = np.empty(len(x_t))
            dy_t[:] = np.nan
            yresid_t = np.empty(len(x_t))
            yresid_t[:] = np.nan
            yadj_t = np.empty(len(x_t))
            yadj_t[:] = np.nan
            self._model._t = self._model.create_t(
                x_t, xmin_t, xmax_t, y_t, dy_t, yresid_t, yadj_t, mask=x_t.mask)
        else:
            self._model._t['X'][where] = x_c
            self._model._t['X'].mask = np.logical_and(
                self._chunk['X'].mask,
                #np.isnan(self._model._t['X']))
                self._model._t['X'].mask)
        y = fun.eval(par, x=np.array(x_c))*np.array(cont_c)
        yresid = y_c-y
        yadj = fun_adj.eval(par_adj, x=np.array(x_c))*cont_c
        self._model._t['Y'][where] = y
        self._model._t['DY'][where] = dy_c
        self._model._t['YRESID'][where] = yresid
        self._model._t['YADJ'][where] = yadj
        self._fun = fun
        self._par = par

        #self.acs.model = self._model


    def N(self, s):
        """ @brief Estimate a column density from an equivalent width 
        @param s A row from a system table
        """

        #\log(N_\mathrm{HI})\simeq 14.096    - 4.6251f   + 18.657f^2 
        #- 46.299f^3 + 53.301f^4 - 23.442f^5

        z = s['Z']
        series = s['SERIES']
        cond_z = self._map['Z'] == z
        

        if series[0:2] == 'Ly' or 'un':

            # Lyman-type series: take Lyman-alpha EW
            if series[0:2] == 'Ly':
                cond_x = self._line.t['X'] == np.max(self._map[cond_z]['X'])
                ew = np.array(self._line.t[cond_x]['EW'])[0]

            # Unknown species: take EW as it were Lyman-alpha's
            else:
                ew = s['EW']
                
            # Parametrization of the curve of growth for Lyman-alpha, b=20
            # To be generalized
            logN_arr = range(12,22)
            lnew_arr = [-2.25, -1.3, -0.7, -0.4, -0.3, -0.2, 0, 0.4, 0.9,\
                         1.4]
            ew_arr = np.exp(lnew_arr)*0.1

            logN = np.min([np.interp(ew, ew_arr, logN_arr), 15])
            s['N'] = 10**logN

    def N_all(self):
        """ @brief Estimate column densities from the equivalent widths 
        """

        for s in self.t:
            self.N(s)

    """
    def extract_residuals(self, kind='abs', diff='max', kappa=3.0):

        spec = dc(self._model)
        print spec.t
        spec.select_extrema(col='YRESID', kind=kind, diff=diff, kappa=kappa)
        print spec._exts_sel
        self._resid_sel = spec._exts_sel
    """        
# To be checked

    def add(self, series, z, **kwargs):
        """ Add a system to a list """

        add = self.create_t(series=series, z=z, **kwargs)
        self._t = vstack([self._t, add])
        #self.create_line()
        #self.group(z)
        #self.chunk(z)
        #self.model(z)
        

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
            for i in dict_series[s['SERIES']]:
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

    def create_t(self, series='unknown',
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
        try:
            t['SERIES'] = Column(series, dtype=object)
        except:
            t.add_column(Column(dtype=object, length=len(t),
                                shape=1, name='SERIES'), index=0)
            for (i, s) in enumerate(t):
                #s['SERIES'] = series[0]
                s['SERIES'] = series
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
            x_temp = self._line.t['X']
            wave = dict_wave[i]
            x = np.append(x, x_temp.to(xunit_def).value)
            ion = np.append(ion, np.full(len(self._line.t), i))
            z = np.append(z, (x_temp/wave).value - 1)
        
        # Redshift table
        self._z = Table()
        self._z['X'] = Column(x, dtype=float, unit=xunit_def)
        self._z['ION'] = Column(ion) 
        self._z['Z'] = Column(z, dtype=float, unit=u.nm/u.nm)

    def extract(self, z, dx=0.0):
        """ Extract a system from a list """
        
        if (hasattr(self, '_group') == False):# or True):
            self.group(z, dx, **kwargs)

        self._t.sort('Z')
        syst_z = self._t['Z']
        syst_idx = np.full(len(self._t), False)
        for g in self._group:
            z = g['Z']
            syst_idx += z==syst_z

        sel = self._t[syst_idx]
        syst_sel = System(
            spec=self._spec, line=self._line, cont=self._cont,
            series=sel['SERIES'],
            z=sel['Z'], N=sel['N'], b=sel['B'], btur=sel['BTUR'],
            dz=sel['DZ'], dN=sel['DN'], db=sel['DB'], dbtur=sel['DBTUR'],
            vary=sel['VARY'], expr=sel['EXPR'])
        return syst_sel, syst_idx
            
    def find(self, series, ztol=1e-4):
        """ Find systems by matching redshifts """
        """ Deprecated """
        
        ion = dict_series[series]
        self.create_z(ion)
        self._z.sort('Z')
        z_arr = self._z['Z']
        ion_arr = self._z['ION']
        match = np.isclose(z_arr[1:], z_arr[:-1], atol=ztol) 
        dec = np.core.defchararray.not_equal(ion_arr[1:], ion_arr[:-1])
        match = np.logical_and(match, dec)
        z_mean = np.mean([z_arr[1:], z_arr[:-1]], axis=0)

        z_sel = z_mean[match]
        self._t = self.create_t(series, z=z_sel)
        self._t['Z'] = z_sel  # To avoid rounding errors
        
        # Table to map rows of self._line.t into rows of self._t
        self._map = Table()
        self._map['X'] = np.append(self._z['X'][1:][match],
                                   self._z['X'][:-1][match])
        self._map['Z'] = np.append(z_sel, z_sel)
        self._map.sort('Z')

        
        
    def merge(self, syst):
        """ Merge two systems """

        self._t = vstack([self._t, syst._t])
        self._t.sort('Z')
        try:
            self._group = vstack([self._group, syst._group])
        except:
            pass
        try:
            self._map = vstack([self._map, syst._map])
        except:
            pass


####
        
    def plot(self, z=None, ax=None, dz=0.008, ions=None):
        """ Plot a system """
        """ Deprecated """

        if (hasattr(self, '_chunk') == False):
            self.chunk(z, **kwargs)

        #t = self._t
        #series = t['SERIES']
        #ions = np.unique([[i for i in s] for s in series])
        rown = 5

        
        """
        if (z == None):
            ions = self._group['ION']
            series = self._group['SERIES']
        else:
            cond = np.logical_and(self._group['Z'] > z-0.002,
                                  self._group['Z'] < z+0.002)
            ions = np.unique(self._group['ION'][np.where(cond)])
            series = self._group['SERIES'][np.where(cond)]
        """
        if (ions is None):
            ions = np.unique([dict_series[i] \
                              for i in self.syst._group['SERIES']])
            ions = ions[np.where(ions != 'unknown')]
        
        #print series
        waves = [dict_wave[i].value for i in ions]
        ions = ions[np.argsort(waves)]
        n = len(ions)
        if ax == None:
            row = min(n,rown)
            col = int(np.ceil(n/rown))
            fig = plt.figure(figsize=(col*6, n*3.5))
            grid = gs(row,col)
            ax = []
            for p in range(n):
                ax.append(fig.add_subplot(grid[p%rown,
                                               int(np.floor(p/rown))]))
            try:
                fig.suptitle(r"$\chi_r^2$ = %3.1f" % self._fit.redchi)
            except:
                pass

        top = []    
        bottom = []
        left = []
        #axt = []
        for p in range(n):
            if p%rown == 0:
                top.append(p)
            if p%rown==rown-1 or p==n-1:
                bottom.append(p)
            if n<rown:
                left.append(p)
            else:
                left = [rown//2]
            #axt.append(ax[p].twiny())
        spec = dc(self._spec.t)
        group = dc(self._group)
        chunk = dc(self._chunk)
        x = spec['X']
        y = spec['Y']
        for c, i in enumerate(ions):
            where_g = group['ION']==i
            zmin = np.min(group['Z'][where_g])-dz
            zmax = np.max(group['Z'][where_g])+dz
            
            x_z = x/dict_wave[i] - 1
            xmin_p = (1+zmin) * dict_wave[i].value
            xmax_p = (1+zmax) * dict_wave[i].value

            #print np.array(group['SERIES']), np.array(series)
            #where_s = [g not in np.array(series) for g in group['SERIES']]
            where_s = [g not in np.array(ions) for g in group['SERIES']]
            #print where_s
            #z = group['Z'][where_g]
            xmin = group['XMIN'][where_g][0]
            xmax = group['XMAX'][where_g][0]

            where_c = np.where(np.logical_and(chunk['X']>=xmin_p,
                                              chunk['X']<=xmax_p))
            
            xc = chunk['X'][where_c]
            yc = chunk['Y'][where_c]
            dyc = chunk['DY'][where_c]
            contc = chunk['CONT'][where_c]
            try:
                modelc = chunk['MODEL'][where_c]
                residc = yc-modelc
            except:
                pass
            xc_z = xc / dict_wave[i] - 1
            
            ax[c].set_xlim(zmin, zmax)
            axt = ax[c].twiny()
            axt.set_xlim(xmin_p, xmax_p)
            maxf = 1.25
            ax[c].set_ylim(-max(contc)*0.2, max(contc)*maxf)
            if c in bottom:
                ax[c].set_xlabel(r"Redshift")
            else:
                ax[c].set_xticks([], [])
            axt.set_xlabel(r"Wavelength [%s]" % x.unit)
            if c in left:
                ax[c].set_ylabel(r"Flux density")#[%s]" % y.unit)
            axt.tick_params(axis="x",direction="in", pad=-20)
            ax[c].text(0.05, 0.5, i, transform=ax[c].transAxes,
                           fontsize=13)
            ax[c].plot(x_z, y, color='C0', linestyle='--')
            #ax[c].plot(x_z, self._conv.y, color='C2', linestyle=':')
            ax[c].plot(xc_z, yc, color='C0')
            ax[c].plot(xc_z, contc, color='C6')
            try:
                ax[c].plot(xc_z, modelc, color='C1')
                ax[c].plot(xc_z, residc, color='C3', linestyle=':')
            except:
                pass
            ax[c].plot(xc_z, dyc, color='C3')
            ax[c].plot(xc_z, -dyc, color='C3')
            #print group
            for g in group[where_g]:
                """
                if c in top:
                    ax[c].text(g['Z'], max(contc)*(maxf+0.18), 
                                   "%3.1f" % np.log10(g['N']), ha='center',
                                   fontsize=9) 
                    ax[c].text(g['Z'], max(contc)*(maxf+0.1), 
                                   "%3.1f" % g['B'], ha='center', fontsize=9) 
                    ax[c].text(g['Z'], max(contc)*(maxf+0.02), 
                                   "%3.1f" % g['BTUR'], ha='center', fontsize=9)
                """
                ax[c].axvline(x=g['Z'], color='C3', alpha=0.5)
                ax[c].axvline(x=g['ZMIN'], color='C3', alpha=0.5,
                              linestyle=':')
                ax[c].axvline(x=g['ZMAX'], color='C3', alpha=0.5,
                              linestyle=':')
            for g in group[where_s]:
                #print where_s
                #print g['X'], g['Z']
                axt.axvline(x=g['X'], color='C3', alpha=0.5, 
                            linestyle='--')
                
        if ax == None:
            grid.tight_layout(fig, rect=[0.01, 0.01, 1, 0.9])
            grid.update(wspace=0.2, hspace=0.0)
            plt.show()

    def save(self, filename):
        #hdu = fits.BinTableHDU.from_columns(
        #    [self._t[i] for i in self._t.columns])
        #hdu.writeto(filename, overwrite=True)
        self._t.write(filename, format='fits', overwrite=True)
            
        
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
        print self_temp._line
        print self_temp._t
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
