from .functions import nfwhm_voigt
from .model import Model
from .vars import *
from astropy import table as at
from astropy import units as au
from astropy import constants as ac
from collections import OrderedDict
from matplotlib import pyplot as plt
import numpy as np

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
        self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))  # Full x array, without NaNs
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
        self._t['pars'] = np.empty(len(self.z), dtype=object)
        self._t['dpars'] = np.empty(len(self.z), dtype=object)
        self._t['chi2r'] = np.empty(len(self.z), dtype=dtype)
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


    def _domain(self, comp, pars, thres=1e-2): #series, z, zmin, zmax):
        """ @brief Define domain for fitting. """

        m = comp.eval(x=self._xs, params=pars)
        c = np.where(m<1-thres)
        xc = np.array(self._xs[c])

        spec = self._spec
        if 'deabs' in spec._t.colnames:
            yc = np.array(spec._t['deabs'][c]/spec._t['cont'][c])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        wc = np.array(spec._t['cont'][c]/spec.dy[c])
        return xc, yc, wc

    def _group(self, comp, pars, thres=1.e-3):
        """ @brief Define group of systems. A group is the set of systems with
        at list one line within the footprint of the system described by the
        model in input.
        """

        # Wavelengths of the lines of all systems
        x_val = [[(1+z.value)*xem_d[t].to(au.nm).value for t in series_d[s]] \
                 for (z,s) in zip(self.z, self.series)]

        x = self._xs
        m = comp.eval(x=x, params=pars)
        group = np.where([np.any(np.interp(xv, x, m)<1-thres) for xv in x_val])
        return group

    def _model_voigt(self, series='Ly_a', z=2.0, N=1e13, b=10, btur=0,
                     resol=35000):

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
        self._t.add_row([series, func, z, dz, zmin, zmax, pars, dpars, chi2r])

        model = Model(series=series, z=z, zmin=zmin, zmax=zmax, pars=pars)
        comp = model._comp
        pars = comp._pars
        return comp, pars

    def _model_voigt_group(self, comp, pars, group):
        for i, s in enumerate(self._t[group[1:]]):
            model = Model(series=s['series'], z=s['z'], zmin=s['zmin'],
                          zmax=s['zmax'], pars=s['pars'], count=i)
            comp *= model._comp
            pars.update(model._comp._pars)
        return comp, pars

    def _update_voigt(self, fit, group=-1):
        """ @brief Update tables after fitting """

        pars = fit.params

        # Spectrum table
        spec = self._spec
        y = spec.y
        cont = spec._t['cont']
        #model = spec._t['model']
        #deabs = spec._t['deabs']
        if 'model' not in spec._t.colnames:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*y.unit
            spec._t['model'][self._s] = cont[self._s]
        spec._t['model'][self._s] = fit.eval(x=self._xs) \
                                    * spec._t['model'][self._s]
        if 'deabs' not in spec._t.colnames:
            spec._t['deabs'] = y
        spec._t['deabs'][self._s] = cont[self._s] + y[self._s] \
                                    - spec._t['model'][self._s]

        # System table
        for i, w in enumerate(group[0]):
            self._t[w]['z'] = pars['lines_voigt_'+str(i)+'_z'].value \
                         #*au.dimensionless_unscaled
            self._t[w]['dz'] = pars['lines_voigt_'+str(i)+'_z'].stderr\
                          #*au.dimensionless_unscaled
            self._t[w]['pars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].value,
                         'b': pars['lines_voigt_'+str(i)+'_b'].value,
                         'btur': pars['lines_voigt_'+str(i)+'_btur'].value,
                         'resol': pars['psf_gauss_'+str(i)+'_resol'].value}
                         #'ampl': pars['adj_gauss_'+str(i)+'_ampl'].value}
            self._t[w]['dpars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].stderr,
                          'b': pars['lines_voigt_'+str(i)+'_b'].stderr,
                          'btur': pars['lines_voigt_'+str(i)+'_btur'].stderr,
                          'resol': pars['psf_gauss_'+str(i)+'_resol'].stderr}
                          #'ampl': pars['adj_gauss_'+str(i)+'_ampl'].stderr}
            self._t[w]['chi2r'] = fit.redchi

    def fit_range(self, series='Ly_a', z_start=2.5, z_end=2.0, z_step=-0.001,
                  thres=1e-3):
        """ @brief Fit the same Voigt model many times in a redshift range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift z_step
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        thres = float(thres)

        for z in np.arange(z_start, z_end, z_step):
            self.fit(series=series, z=z, thres=thres)
            print(prefix, "I've fitted system at redshift %2.4f…" % z, end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f with a step of %2.4f." % (z_start, z_end, z_step))

        return 0


    def fit(self, series='Ly_a', z=2.0, N=1e13, b=10, btur=0, resol=35000,
            ampl=0.0, thres=1e-3):
        """ @brief Create and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Redshift
        @param N Column density
        @param b Doppler broadening
        @param btur Turbulence broadening
        @param resol Resolution
        @param ampl Amplitude of the continuum adjustment
        @param thres Threshold for grouping
        @return 0
        """

        z = float(z)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)
        ampl = float(ampl)
        thres = float(thres)

        comp, pars = self._model_voigt(series, z, N, b, btur, resol)
        group = self._group(comp, pars, thres)
        comp, pars = self._model_voigt_group(comp, pars, group)
        xc, yc, wc = self._domain(comp, pars, thres)
        fit = comp.fit(yc, pars, x=xc, weights=wc)
        self._update_voigt(fit, group)

        return 0
