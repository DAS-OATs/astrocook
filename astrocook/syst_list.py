from .model import Model
from .vars import *
from astropy import table as at
from astropy import units as au
from collections import OrderedDict
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


    def _domain(self, series, zmin, zmax):
        """ @brief Define domain for fitting. """

        c = np.array([], dtype=int)
        for t in series_d[series]:
            c = np.append(c, np.where(
                    np.logical_and(self._xs > (1+zmin)*xem_d[t].value,
                                   self._xs < (1+zmax)*xem_d[t].value))[0])
        xc = np.array(self._xs[c])  # Intervals used for fitting

        spec = self._spec
        if 'deabs' in spec._t.colnames:
            yc = np.array(spec._t['deabs'][c]/spec._t['cont'][c])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        wc = np.array(spec._t['cont'][c]/spec.dy[c])
        return xc, yc, wc

    def _find_group(self, z_sel, x_tol=0.000001):
        """ @brief Find group of systems around a given redshift. A group is the
        set of all systems that have at least one line closer than the tolerance
        to at least one line of the system at the given redshift.
        """

        # Wavelengths of the lines of the system at the given redshift
        sel = self.z == z_sel
        x_sel = [(1+z_sel)*xem_d[t].to(au.nm).value \
                 for t in series_d[self.series[sel][0]]]

        # Wavelengths of the lines of all systems
        x_val = [[(1+z.value)*xem_d[t].to(au.nm).value for t in series_d[s]] \
                 for (z,s) in zip(self.z, self.series)]

        # Absolute difference between the first two
        x_abs = [[[np.absolute(xs-xv) for xs in x_sel] for xv in x_val[i]] \
                 for i in range(len(x_val))]

        # Array to select systems with absolute difference below tolerance
        self._group = [np.any(np.array(xa)<x_tol) for xa in x_abs]


    def _update_voigt(self, fit, group=-1):
        """ @brief Update tables after fitting """

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
        #bv = fit.best_values
        pars = fit.params
        where = np.where(np.array(self._group) == 1)[0]
        for i, w in enumerate(where):
            self._t[w]['z'] = pars['lines_voigt_'+str(i)+'_z'].value \
                         #*au.dimensionless_unscaled
            self._t[w]['dz'] = pars['lines_voigt_'+str(i)+'_z'].stderr\
                          #*au.dimensionless_unscaled
            self._t[w]['pars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].value,
                         'b': pars['lines_voigt_'+str(i)+'_b'].value,
                         'btur': pars['lines_voigt_'+str(i)+'_btur'].value,
                         'resol': pars['psf_gauss_'+str(i)+'_resol'].value}
            self._t[w]['dpars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].stderr,
                          'b': pars['lines_voigt_'+str(i)+'_b'].stderr,
                          'btur': pars['lines_voigt_'+str(i)+'_btur'].stderr,
                          'resol': pars['psf_gauss_'+str(i)+'_resol'].stderr}
            self._t[w]['chi2r'] = fit.redchi

    def fit_many(self, series='Ly_a', z_start=2.5, z_end=2.0, z_step=-0.01):
        """ @brief Fit the same Voigt model many times in a redshift range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift z_step
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)

        for z in np.arange(z_start, z_end, z_step):
            self.fit_single(series=series, z=z)
            print(prefix, "I've fitted system at redshift %2.4f…" % z, end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f" % (z_start, z_end))

        return 0


    def fit_single(self, series='Ly_a', z=2.0, z_tol=0.002, N=1e13, b=10,
                   btur=0, resol=35000, cont_ampl=0.1):
        """ @brief Create and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Redshift
        @param z_tol Redshift tolerance
        @param N Column density
        @param b Doppler broadening
        @param btur Turbulence broadening
        @param resol Resolution
        @return 0
        """

        z = float(z)
        z_tol = float(z_tol)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)
        cont_ampl = float(cont_ampl)

        # Create model
        func = 'voigt'
        dz = 0.0
        zmin = z-z_tol
        zmax = z+z_tol
        pars = {'N': N, 'b': b, 'b_tur': btur, 'resol': resol}
        dpars = {'N': None, 'b': None, 'b_tur': None, 'resol': None}
        chi2r = None
        self._t.add_row([series, func, z, dz, zmin, zmax, pars, dpars, chi2r])

        self._find_group(z)
        """
        model = Model(series=series, z=z, zmin=zmin, zmax=zmax, pars=pars)
        comp = model._comp
        comp._pars.pretty_print()
        """
        for i, s in enumerate(self._t[self._group]):
            model = Model(series=s['series'], z=s['z'], zmin=s['zmin'],
                          zmax=s['zmax'], pars=s['pars'], count=i)
            if i == 0:
                comp = model._comp
                #comp._pars.pretty_print()
                pars = model._comp._pars
            else:
                comp *= model._comp
                #comp._pars.pretty_print()
                pars.update(model._comp._pars)
        #"""

        # Fit model
        xc, yc, wc = self._domain(series, zmin, zmax)
        #fit = comp.fit(yc, comp._pars, x=xc, weights=wc)
        fit = comp.fit(yc, pars, x=xc, weights=wc)
        self._update_voigt(fit)

        return 0
