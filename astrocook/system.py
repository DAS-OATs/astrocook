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
                 z=None,
                 #N=N_def,
                 N=None,
                 b=b_def, btur=btur_def,  
                 dz=None, dN=None, db=None, dbtur=None,
                 vary=[True,True,True,False], expr=[None,None,None,None],
                 chi2r=None,
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

        # Continuum
        if (cont != None):
            self._cont = dc(cont)

            
        # "is not" works also with arrays
        if (z is not None):
            if (btur is None):
                btur = 0.0 * bunit_def
            self._t = self.create_t(series, z, N, b, btur,
                                    dz, dN, db, dbtur, vary, expr, chi2r,
                                    Nunit, bunit, dtype)
            
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
        try:
            self._map = acs.line._map  # When loading from a current session
        except:
            pass  # When loading from a saved session
        
        
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
        #print len(where)

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
        if s != None:
            self.group(s)
            self.chunk()
        
        out = self._model.apply_mask(
            ['X', 'X'], [range(len(self._model.t)), self._chunk_rows],
            [True, False])
        out.t['Y'] = out.t['YRESID']
        out.t['Y'][out.t['X'].mask] = np.nan
        out.t.remove_columns(['YRESID', 'YADJ'])
        return out

        
    def fit(self, s=None, z=None, done=False, **kwargs):
        """ @brief Fit a model on a group of lines 
        
        @param s A row from a system table
        """

        if z != None:
            s = self._t[self._t['Z']==z]
        
        if s == None:
            try:
                s = self._syst_sel
            except:
                s = self.acs._syst_sel
        
        self.model(s, **kwargs)
        #print self._group['DONE']

        where = self._chunk['X'].mask == False
        x_c = self._chunk['X'][where]
        y_c = self._chunk['Y'][where]
        dy_c = self._chunk['DY'][where]
        cont_c = self._chunk['CONT'][where]

        fun = self._fun
        par = self._par
        #par.pretty_print()
        fit = fun.fit(np.array(y_c/cont_c), par, x=np.array(x_c),
                      weights=np.array(cont_c/dy_c))
        par = fit.params
        #par.pretty_print()
        y = fit.eval(par, x=x_c) * cont_c
        yresid = y_c-y
        yadj = fit.eval_components()['adj_'] * cont_c
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
            pref = l['PREF']
            #print pref, l['Z'], par[pref+'_z'].value
            if par[pref+'_z'].value > l['ZMIN'] and \
               par[pref+'_z'].value < l['ZMAX'] or 1==1:
                try:
                    l['Z'] = par[pref+'_z'].value
                    l['N'] = par[pref+'_N'].value
                    l['B'] = par[pref+'_b'].value
                    l['BTUR'] = par[pref+'_btur'].value
                    l['DZ'] = par[pref+'_z'].stderr
                    l['DN'] = par[pref+'_N'].stderr
                    l['DB'] = par[pref+'_b'].stderr
                    l['DBTUR'] = par[pref+'_btur'].stderr
                    l['CHI2R'] = fit.redchi 
                    l['DONE'] = done 
                    l['X'] = (1.+l['Z'])*dict_wave[l['ION']].value
                except:
                    print "Group not updated."
            else:
                l_del.append(i)
        #print self._group['DONE']
        
        # Save fitted lines in the group
        try:
            ok
        except:
            cond = self._group['Z'] > 0 
            try:
                new_t = unique(self._group[self._t.colnames][cond], keys='Z')
                #print new_t
                #print self._t[self._group_t]
                self._t[self._group_t] = new_t#[new_t['Z']>0]
            except:
                new_t = unique(self._group[self._t.colnames][cond], keys='N')
                #print new_t
                #print self._t[self._group_t]
                self._t[self._group_t] = new_t#[new_t['Z']>0]
            new_map = self._group[self._map.colnames][cond]
            new_line = self._group[self._line.t.colnames]
            self._map[self._group_map] = new_map
            self._line._t[self._group_line] = new_line
            self._line._t.sort('X')
        #except:
        #    print "Lists not updated."

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

        self._line.t.sort('X')

        # Join systems and lines
        join_t = join(join(self._t, self._map), self._line.t)
        #print self._map[np.logical_and(self._map['X']>370.862, self._map['X']<372.905)]
        #print join_t[np.logical_and(join_t['Z']>2.057, join_t['Z']<2.060)]
        #print self._t[np.logical_and(self._t['Z']>2.057, self._t['Z']<2.060)]
        #print self._map[np.logical_and(self._map['Z']>2.057, self._map['Z']<2.060)]
        
        # Select the system redshift 
        join_z = join_t['Z']
        cond_z = s['Z']==join_z
        #print s['Z']
        # Select other systems close in redshift
        diff = 1
        while diff > 0:
            join_xmin = join_t['XMIN']
            join_xmax = join_t['XMAX']
            cond_x = np.full(len(join_t), False)
            for j in join_t[cond_z]:
                xmin = j['XMIN']
                xmax = j['XMAX']
                cond_x += np.logical_and(join_xmax>xmin, join_xmin<xmax)

            diff = np.abs(np.sum(cond_x)-np.sum(cond_z))

            # ...then select the whole systems those lines belong to
            cond_z = np.full(len(join_t), False)
            for j in join_t[np.where(cond_x)[0]]:
                z = j['Z']
                cond_z += z==join_z
            
        group = join_t[cond_z]
        #print group
        
        # Find wavelength duplicates (it may happen that a line may be
        # associated to two systems as the same ion) and remove them from group
        # and map
        group.sort('X')
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
        if len(zlist)>0:
            diff1d = np.append(zlist[0], np.ediff1d(zlist))  
            where = np.where(diff1d == 0)[0]

            # Associate each duplicate row to its companion, to link parameters 
            for (l, w) in enumerate(where):

                # Change the expression used in fit
                p = group['PREF'][w-1]
                group['VARY'][w] = [False, False, False, False]
                group['EXPR'][w] = [p+'_z', p+'_N', p+'_b', p+'_btur']

        self._group = group
        #print self._group
        #print self._line.t[315:325]
        
        self._t.sort('Z')
        self._map.sort('Z')
        self._group_t = np.in1d(self._t['Z'], group['Z'])
        self._group_map = np.in1d(self._map['Z'], group['Z'])
        self._group_line = np.in1d(self._line.t['X'], group['X'])

    def line_merge(self, series='Ly_ab', keep='complete'):
        """ @brief Merge new lines into a list of systems

        @param series Label of the series of transitions
        @param keep 'all': all line redshifts are defined as systems;
                    'match': only matching redshifts are defined as systems;
                    'complete': matching redshifts are defined as systems, and
                    then remaining line redshifts are added to the system list
        """ 

        # Create a new list of systems from the line list
        temp = self.line_new(series, keep)
        temp._map.sort('X')

        # Find only the new lines in the list
        try:
            new = self.acs.line._new

            # Update the map
            null1, null2, int_map = np.intersect1d(new['X'], temp._map['X'],
                                                   return_indices=True)
            #new_map = temp._map[new]
            self._map = vstack([self._map, temp._map[int_map]])
            #self._map = unique(self._map, keys='Z')
            self._map.sort('Z')
            #print "temp_t", np.array(temp._t[np.logical_and(temp._t['Z']>2.32131, temp._t['Z']<2.32133)])
            
            # Update the system list adding only the new systems
            null1, null2, int_t = np.intersect1d(
                temp._map['Z'][int_map], temp._t['Z'], return_indices=True)
            #print temp._t[int_t]
            self._t = vstack([self._t, temp._t[int_t]])

            self._map = unique(self._map, keys='X')
            self._map.sort('Z')

            self._t = unique(self._t, keys='Z')
            self._t.sort('Z')
            #return self._t
        except:
            pass

        #print "self_t", np.array(self._t[np.logical_and(self._t['Z']>2.32131, self._t['Z']<2.32133)])
        
    def line_new(self, series='Ly_ab', keep='complete'):
        """ @brief Use matching redshift from a list of lines to create a list
        of systems

        @param series Label of the series of transitions
        @param keep 'all': all line redshifts are defined as systems;
                    'match': only matching redshifts are defined as systems;
                    'complete': matching redshifts are defined as systems, and
                    then remaining line redshifts are added to the system list
        """ 

        
        line = self.acs.line
        line_map = Table()
        if keep == 'all':
            z = np.array(line._z['Z'])
            line_map['X'] = line._z['X']
            line_map['Z'] = line._z['Z']
        if keep == 'match' or keep == 'complete':
            z = line._z_match
            #print "z", z[np.logical_and(z>2.05875, z<2.058755)]
            line_map['X'] = np.append(line._z['X'][1:][line._w_match],
                                      line._z['X'][:-1][line._w_match])
            line_map['Z'] = np.append(line._z_match, line._z_match)
        out = System(acs=self.acs, series=series, z=z)
        line_map_add = Table()
        
        # If series is Ly_ab(...), discarded lines are added as Ly_a's
        if keep == 'complete':
            disc = np.logical_and(~np.in1d(line._z['X'], line_map['X']),
                                  line._z['ION'] == dict_series[series][-1])
            z_add = line._z['Z'][disc]
            #print "z_add", np.array(z_add[np.logical_and(z_add>2.32131, z_add<2.32133)])
            line_map_add['X'] = line._z['X'][disc]
            line_map_add['Z'] = line._z['Z'][disc]
            line_map = vstack([line_map, line_map_add])
            out_add = System(acs=self.acs, series=dict_series[series][-1],
                             z=z_add)
            out._t = vstack([out._t, out_add._t])
            out._t.sort('Z')
        line_map.sort('Z')
        out._map = line_map 
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
        
        mod = Model()

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
            prof_vary = l['VARY']
            prof_expr = l['EXPR']
            prof_pref = l['PREF']
            getattr(mod, prof)(ion, wave, value=prof_value,
                               min=[l['ZMIN'], 1e10, 0.0, None],
                               max=[l['ZMAX'], 1e23, 200.0, None],
                               vary=prof_vary,
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
        if len(rem) == 0:
            print "Can't model this system: chunk empty."
            
        for i, (c_min, c_max, c_mean) in enumerate(zip(imin, imax, imean)):
            #print i, c_min, c_max, c_mean, rem['X'][c_mean]
            #print psf_resol_value
            #if psf_resol_value < 0.0 or 1==1:
            #if len(rem) > 0:
            psf_resol_value = rem['RESOL'][c_mean]
            #print psf_resol_value
            psf_value = [c_min, c_max, rem['X'][c_mean], psf_resol_value]
            psf_vary = [False, False, False, psf_resol_vary]
            psf_min = [None, None, None, psf_resol_min]
            psf_max = [None, None, None, psf_resol_max]
            psf_expr = [None, None, None, psf_resol_expr]
            psf_pref = 'psf_'+str(i) 
            getattr(mod, psf)(psf_value, psf_vary, psf_min, psf_max,
                              psf_expr, psf_pref)
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

        z = s['Z']
        series = s['SERIES']
        cond_z = self._map['Z'] == z
        if (series[0:2] == 'Ly' or 'un') and np.sum(cond_z)>0:

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
        else:
            s['N'] = 1e14

    def N_all(self):
        """ @brief Estimate column densities from the equivalent widths 
        """

        #print len(self.t)
        for s in self.t:
            if np.isnan(s['N']):
                self.N(s)

# To be checked

    #"""
    def add(self, series, z, **kwargs):
        
        add = self.create_t(series=series, z=z, **kwargs)
        self._t = vstack([self._t, add])
        #self.create_line()
        #self.group(z)
        #self.chunk(z)
        #self.model(z)
        

    def create_line(self, xmin=None, xmax=None, sigma=0.07):
        
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
                 z=z_def,
                 #N=N_def,
                 N=None,
                 b=b_def, btur=btur_def,
                 dz=None, dN=None, db=None, dbtur=None,
                 vary=[True, True, True, False], expr=[None,None,None,None],
                 chi2r=None,
                 Nunit=Nunit_def, bunit=bunit_def, dtype=float):
        
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
        chi2r = np.array(chi2r, ndmin=1)
        done = np.array(False, ndmin=1)
        zunit = u.nm/u.nm
        #print "create_t z", z[np.logical_and(z>2.05875, z<2.058755)]
        t = Table()
        t['Z'] = Column(z, dtype=dtype, unit=zunit)
        #print "create_t tz", np.array(t['Z'][np.logical_and(t['Z']>2.05875, t['Z']<2.058755)])
        t['N'] = Column(N, dtype=dtype, unit=Nunit)
        t['B'] = Column(b, dtype=dtype, unit=bunit)
        t['BTUR'] = Column(btur, dtype=dtype, unit=bunit)
        t['DZ'] = Column(dz, dtype=dtype, unit=zunit)
        t['DN'] = Column(dN, dtype=dtype, unit=Nunit)
        t['DB'] = Column(db, dtype=dtype, unit=bunit)
        t['DBTUR'] = Column(dbtur, dtype=dtype, unit=bunit)
        t['VARY'] = Column(vary, dtype=object)
        t['EXPR'] = Column(expr, dtype=object)  # dtype to avoid truncation
        t['CHI2R'] = Column(chi2r, dtype=dtype)
        t['DONE'] = Column(done, dtype=bool)
        
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
