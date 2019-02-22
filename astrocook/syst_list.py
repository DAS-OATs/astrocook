from .model import Model, ModelLines, ModelPSF
from .spectrum import Spectrum
from .vars import *
from .functions import convolve, lines_voigt, psf_gauss
from astropy import table as at
from astropy import units as au
from astropy import constants as ac
from collections import OrderedDict
from copy import deepcopy as dc
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import argrelmax

prefix = "System list:"

class SystList(object):
    """ Class for system lists

    A SystList is a list of absorption systems with methods for handling
    spectral lines. """

    def __init__(self,
                 sess,
                 series=[],
                 func=[],
                 z=[],
                 dz=[],
                 zmin=[],
                 zmax=[],
                 dtype=float):
        self._spec = sess.spec
        self._lines = sess.lines
        self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))  # Full x array, without NaNs
        self._ys = np.ones(len(self._xs))
        self._s = self._spec._where_safe

        t = at.Table()
        zunit = au.dimensionless_unscaled
        t['series'] = at.Column(np.array(func, ndmin=1), dtype='S100')
        t['func'] = at.Column(np.array(func, ndmin=1), dtype='S5')
        t['z'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        t['dz'] = at.Column(np.array(dz, ndmin=1), dtype=dtype, unit=zunit)
        t['zmin'] = at.Column(np.array(zmin, ndmin=1), dtype=dtype, unit=zunit)
        t['zmax'] = at.Column(np.array(zmax, ndmin=1), dtype=dtype, unit=zunit)
        self._t = t
        self._t['comp'] = np.empty(len(self.z), dtype=object)
        self._t['pars'] = np.empty(len(self.z), dtype=object)
        self._t['dpars'] = np.empty(len(self.z), dtype=object)
        self._t['chi2r'] = np.empty(len(self.z), dtype=dtype)
        self._t['ys'] = np.empty(len(self.z), dtype=object)
        self._dtype = dtype


    @property
    def t(self):
        return self._t

    @property
    def series(self):
        return self._t['series']

    @property
    def func(self):
        return self._t['func']

    @property
    def z(self):
        return au.Quantity(self._t['z'])

    @property
    def zmin(self):
        return au.Quantity(self._t['zmin'])

    @property
    def zmax(self):
        return au.Quantity(self._t['zmax'])

    @series.setter
    def func(self, val):
        self._t['series'] = np.array(val, dtype='S100')

    @func.setter
    def func(self, val):
        self._t['func'] = np.array(val, dtype='S5')

    @z.setter
    def z(self, val, dtype=float):
        self._t['z'] = np.array(val, dtype=dtype)
        self._t['z'].unit = val.unit

    @zmin.setter
    def zmin(self, val, dtype=float):
        self._t['zmin'] = np.array(val, dtype=dtype)
        self._t['zmin'].unit = val.unit

    @zmax.setter
    def zmax(self, val, dtype=float):
        self._t['zmax'] = np.array(val, dtype=dtype)
        self._t['zmax'].unit = val.unit


    #def _domain(self, comp, pars, thres=1e-3): #series, z, zmin, zmax):
    def _domain(self, ys, thres=1e-3): #series, z, zmin, zmax):
        """ @brief Define domain for fitting. """

        spec = self._spec
        #m = comp.eval(x=self._xs, params=pars)
        #c = np.where(m<1-thres)
        c = np.where(ys<1-thres)
        xc = np.array(self._xs[c])

        if 'deabs' in spec._t.colnames:# and 1 == 0:
            yc = np.array(spec._t['deabs'][c]/spec._t['cont'][c])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        wc = np.array(spec._t['cont'][c]/spec.dy[c])
        return xc, yc, wc

    def _fit_voigt(self, comp, pars, xc, yc, wc):

        fit = comp.fit(yc, pars, x=xc, weights=wc)
        comp = fit
        pars = fit.params
        self._ys = comp.eval(x=self._xs, params=pars)
        return comp, pars

    def _group_voigt(self, comp, pars, ys, thres=1.e-3):
        """ @brief Define group of systems. A group is the set of systems with
        at list one line within the footprint of the system described by the
        model in input.
        """

        group = [len(self._t)-1]
        c = 0
        for i, s in enumerate(self._t[:-1]):
            yn = s['ys']
            if np.amin(np.maximum(ys, yn)) < 1-thres:
                group.append(i)
                #m = Model(series=s['series'], z=s['z'], zmin=s['zmin'],
                #          zmax=s['zmax'], pars=s['pars'], count=c+1)
                m = ModelLines(lines_voigt, s['series'], s['z'], s['zmin'],
                               s['zmax'], c+1, **s['pars'])
                #comp *= m._comp
                #pars.update(m._comp._pars)
                comp *= m
                pars.update(m._pars)

                c += 1

        group = np.array(group)
        ys = comp.eval(x=self._xs, params=pars)
        return comp, pars, ys, group

    def _single_voigt(self, series='Ly_a', z=2.0, N=1e13, b=10, btur=0,
                      resol=35000, add=True):

        # Create model
        func = 'voigt'
        dz = 0.0
        z_tol = 1e-3
        zmin = z-z_tol
        zmax = z+z_tol
        pars = {'N': N, 'b': b, 'b_tur': btur, 'resol': resol}#, 'ampl': ampl}
        dpars = {'N': None, 'b': None, 'b_tur': None, 'resol': None}#,
#                 'ampl': None}
        chi2r = None

        #m = Model(series=series, z=z, zmin=zmin, zmax=zmax, pars=pars)
        m = ModelLines(lines_voigt, series, z, zmin, zmax, 0, **pars)
        #comp = m._comp
        comp = m
        pars = comp._pars
        ys = comp.eval(x=self._xs, params=pars)

        if add:
            self._t.add_row([series, func, z, dz, zmin, zmax, comp, pars, dpars,
                             chi2r, ys])

        return comp, pars, ys

    def _update_spec(self, fit=None):#, fit):
        """ @brief Update spectrum after fitting """

        spec = self._spec
        y = spec.y
        if 'model' not in spec._t.colnames:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*y.unit
        if 'deabs' not in spec._t.colnames:
            spec._t['deabs'] = y

        s = self._s
        cont = spec._t['cont']
        model = spec._t['model']
        deabs = spec._t['deabs']

        model[s] = cont[s]
        for i, r in enumerate(self._t):
            #m = Model(series=r['series'], z=r['z'], zmin=r['zmin'],
            #          zmax=r['zmax'], pars=r['pars'])
            print(r['pars'])
            m = ModelLines(lines_voigt, r['series'], r['z'], r['zmin'],
                           r['zmax'], i, **r['pars'])
            #model[s] = m._comp.eval(x=self._xs, params=m._comp._pars) * model[s]
            model[s] = m.eval(x=self._xs, params=m._pars) * model[s]
        deabs[s] = cont[s] + y[s] - model[s]


    def _update_systs(self, comp, pars, ys, group=None):
        """ @brief Update system list after fitting """

        #pars = fit.params
        if group is None:
            group = [len(self._t)-1]
        for i, w in enumerate(group):
            self._t[w]['comp'] = comp
            self._t[w]['z'] = pars['lines_voigt_'+str(i)+'_z'].value \
                         #*au.dimensionless_unscaled
            self._t[w]['dz'] = pars['lines_voigt_'+str(i)+'_z'].stderr\
                          #*au.dimensionless_unscaled
            #"""
            self._t[w]['pars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].value,
                         'b': pars['lines_voigt_'+str(i)+'_b'].value,
                         'btur': pars['lines_voigt_'+str(i)+'_btur'].value}
                         #'resol': pars['psf_gauss_'+str(i)+'_resol'].value}
                         #'ampl': pars['adj_gauss_'+str(i)+'_ampl'].value}
            self._t[w]['dpars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].stderr,
                          'b': pars['lines_voigt_'+str(i)+'_b'].stderr,
                          'btur': pars['lines_voigt_'+str(i)+'_btur'].stderr}
                          #'resol': pars['psf_gauss_'+str(i)+'_resol'].stderr}
                          #'ampl': pars['adj_gauss_'+str(i)+'_ampl'].stderr}
            #"""
            self._t[w]['chi2r'] = comp.redchi

    def fit_from_lines(self, series='Ly_a', z_start=2.5, z_end=2.0, N=1e14,
                       b=10, thres=1e-3):
        """ @brief Fit Voigt models to a line list, given a redshift range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param N Guess column density
        @param b Guess doppler broadening
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        N = float(N)
        b = float(b)
        thres = float(thres)

        z_lines = [[(x.to(au.nm)/xem_d[t].to(au.nm))-1. \
                    for t in series_d[series]] for x in self._lines.x]
        z_lines = np.ravel(z_lines)
        if z_end < z_start:
            z_range = z_lines[np.logical_and(z_lines<z_start, z_lines>z_end)]\
                          [::-1]
        else:
            z_range = z_lines[np.logical_and(z_lines>z_start, z_lines<z_end)]

        #for z in np.arange(z_start, z_end, z_step):
        for z in z_range:
            self.fit(series=series, z=z, N=N, b=b, group_thres=thres,
                     update=False)
            print(prefix, "I've fitted a system at redshift %2.4f…" % z,
                  end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f." % (z_start, z_end))

        self._update_spec()

        return 0

    def fit_from_deabs(self, series='Ly_a', z_start=2.5, z_end=2.0, N=1e12,
                       b=10, thres=1e-3):
        """ @brief Fit Voigt models to spectrum residuals, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param N Guess column density
        @param b Guess doppler broadening
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        N = float(N)
        b = float(b)
        thres = float(thres)

        spec_deabs = dc(self._spec)
        spec_deabs.convolve_gauss(input_col='deabs')
        lines_deabs = spec_deabs.find_peaks()

        z_lines = [[(x.to(au.nm)/xem_d[t].to(au.nm))-1. \
                    for t in series_d[series]] for x in lines_deabs.x]
        z_lines = np.ravel(z_lines)
        if z_end < z_start:
            z_range = z_lines[np.logical_and(z_lines<z_start, z_lines>z_end)]\
                          [::-1]
        else:
            z_range = z_lines[np.logical_and(z_lines>z_start, z_lines<z_end)]

        #for z in np.arange(z_start, z_end, z_step):
        for z in z_range:
            self.fit(series=series, z=z, N=N, b=b, group_thres=thres,
                     update=False)
            print(prefix, "I've fitted a system at redshift %2.4f…" % z,
                  end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f." % (z_start, z_end))

        self._update_spec()

        return 0

    def fit_range(self, series='Ly_a', z_start=2.5, z_end=2.0, z_step=1e-3,
                  N=1e13, b=10, thres=1e-3):
        """ @brief Fit Voigt models at constant steps in redshift range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift step
        @param N Guess column density
        @param b Guess doppler broadening
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        N = float(N)
        b = float(b)
        thres = float(thres)

        for z in np.arange(z_start, z_end, z_step):
            self.fit(series=series, z=z, N=N, b=b, group_thres=thres,
                     update=False)
            print(prefix, "I've fitted a system at redshift %2.4f…" % z,
                  end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f." % (z_start, z_end))

        self._update_spec()

        return 0

    def fit(self, series='Ly_a', z=2.0, N=1e13, b=10, btur=0, resol=35000,
            ampl=0.0, group_thres=1e-3, domain_thres=1e-3, update=True):
        """ @brief Create and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param btur Guess turbulence broadening
        @param resol Resolution
        @param ampl Amplitude of the continuum adjustment
        @param group_thres Threshold for grouping
        @param domain_thres Threshold for fitting domain
        @param update Flag to update the spectrum
        @return 0
        """

        z = float(z)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)
        ampl = float(ampl)
        group_thres = float(group_thres)
        domain_thres = float(domain_thres)

        comp, pars, ys = self._single_voigt(series, z, N, b, btur, resol)
        comp, pars, ys, group = self._group_voigt(comp, pars, ys, group_thres)
        #xc, yc, wc = self._domain(comp, pars, domain_thres)
        xc, yc, wc = self._domain(ys, domain_thres)
        #comp, pars = self._conv_psf(comp, pars)
        comp, pars = self._fit_voigt(comp, pars, xc, yc, wc)
        self._update_systs(comp, pars, group)
        if update:
            self._update_spec(comp)

        return 0
