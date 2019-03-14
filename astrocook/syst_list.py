#from .model import Model, ModelLines, ModelPSF
from .model_list import ModelList
from .syst_model import SystModel
#from .spectrum import Spectrum
from .vars import *
from .functions import convolve, lines_voigt, psf_gauss, running_mean
from astropy import table as at
from astropy import units as au
#from astropy import constants as ac
#from astropy.stats import sigma_clip
#from collections import OrderedDict
from copy import deepcopy as dc
#from lmfit import CompositeModel as LMComposite
from matplotlib import pyplot as plt
import numpy as np
#from scipy.signal import argrelmax

prefix = "System list:"

class SystList(object):
    """ Class for system lists

    A SystList is a list of absorption systems with methods for handling
    spectral lines. """

    def __init__(self,
#                 sess,
                 spec=None,
                 sess=None,
                 series=[],
                 func=[],
                 z=[],
                 logN=[],
                 b=[],
                 dtype=float):
        if sess != None:
            self._spec = sess.spec
            self._lines = sess.lines
        else:
            self._spec = spec #sess.spec
            self._lines = spec._lines #sess.lines
        self._mods = ModelList()#sess)
        self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))  # Full x array, without NaNs
        self._ys = np.ones(len(self._xs))
        self._s = self._spec._where_safe

        t = at.Table()
        zunit = au.dimensionless_unscaled
        logNunit = au.dimensionless_unscaled
        bunit = au.km/au.s
        t['func'] = at.Column(np.array(func, ndmin=1), dtype='S5')
        t['series'] = at.Column(np.array(func, ndmin=1), dtype='S100')
        t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        t['z'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        t['logN'] = at.Column(np.array(logN, ndmin=1), dtype=dtype,
                              unit=logNunit)
        t['b'] = at.Column(np.array(b, ndmin=1), dtype=dtype, unit=bunit)
        self._t = t
        #self._t['mod'] = np.empty(len(self.z), dtype=object)
        self._t['chi2r'] = np.empty(len(self.z), dtype=dtype)
        self._t['count'] = np.empty(len(self.z), dtype=int)
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

    @series.setter
    def series(self, val):
        self._t['series'] = np.array(val, dtype='S100')

    @func.setter
    def func(self, val):
        self._t['func'] = np.array(val, dtype='S5')

    @z.setter
    def z(self, val, dtype=float):
        self._t['z'] = np.array(val, dtype=dtype)
        self._t['z'].unit = val.unit

    def _add(self, series='Ly_a', z=2.0, logN=13, b=10, resol=70000):
        """ @brief Add a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @return 0
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        vars = {'z': z, 'logN': logN, 'b': b, 'resol': resol}
        count = len(self._t)

        self._mod = SystModel(self, series, vars, z0=z)
        self._t.add_row(['voigt_func', series, z, z, logN, b, None, count])

        return 0

    def _add_fit(self, series='Ly_a', z=2.0, logN=14, b=10, resol=70000,
                 chi2r_thres=None, verb=False):
        """ @brief Add a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @return 0
        """

        z_range = z if np.size(z) > 1 else [z]
        for i, z in enumerate(z_range):
            if verb:
                print(prefix, "I'm fitting a %s system at redshift %2.4f "\
                      "(%i/%i)…" % (series, z, i+1, len(z_range)), end='\r')
            self._add(series, z, logN, b, resol)
            self._fit(chi2r_thres)
        if verb:
            print(prefix, "I've fitted %i %s systems between redshift %2.4f "\
                  "and %2.4f." % (len(z_range), series, z_range[0],
                                  z_range[-1]))
        return 0

    def _append(self, frame):
        vstack = at.vstack([self._t, frame._t])
        self._t = at.unique(vstack, keys=['z0'])
        return 0

    def _fit(self, chi2r_thres=None):
        """ @brief Fit a Voigt model for a system.
        @param mod Model
        @return 0
        """

        if chi2r_thres == None:
            chi2r_thres = np.inf

        self._mod.fit()
        self._mod._pars.pretty_print()
        for row in range(len(self._t)):
            try:
                if self._mod._chi2r < chi2r_thres:
                    count = self._t[row]['count']
                    pref = 'lines_voigt_'+str(count)
                    self._t[row]['z'] = self._mod._pars[pref+'_z']
                    self._t[row]['logN'] = self._mod._pars[pref+'_logN']
                    self._t[row]['b'] = self._mod._pars[pref+'_b']
                    self._t[row]['chi2r'] = self._mod._chi2r
                else:
                    where = np.where(self._mods._t['z0'] == self._t[row]['z0'])
                    self._mods._t.remove_rows(where)
                    self._t.remove_row(row)
            except:
                pass

        #print(self._t)
        return 0

    def _test(self, spec, xm, ym, ym_0, ym_1, ym_2, col='y', chi2_fact=0.8):
        """ @brief Test a Voigt model for a system.
        """

        ys = np.interp(xm, spec.x.to(au.nm), spec._t[col]/spec._t['cont'])
        dys = np.interp(xm, spec.x.to(au.nm), spec.dy/spec._t['cont'])
        chi2 = np.sum(((ys-ym)/dys)**2)
        chi2_0 = np.sum(((ys-ym_0)/dys)**2)
        chi2_1 = np.sum(((ys-ym_1)/dys)**2)
        chi2_2 = np.sum(((ys-ym_2)/dys)**2)
        if chi2 < chi2_fact*np.min([chi2_0, chi2_1, chi2_2]):
            return True, chi2, chi2_0
        else:
            return False, chi2, chi2_0

    def _test_fit(self, spec, series='Ly_a', z=2.0, logN=14, b=10, resol=75000,
                  col='y', chi2_fact=1.0, chi2r_thres=2.0, verb=False):

        z_range = z if np.size(z) > 1 else [z]

        spec_temp = spec
        vars = {'z': 0.0, 'logN': logN, 'b': b, 'resol': resol}

        z_min = np.min([(np.min(self._spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_max = np.max([(np.max(self._spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])

        z_range = z_range[np.where(np.logical_and(z_range > z_min,
                                                  z_range < z_max))]
        self._spec._shift_rf(0.5*(z_min+z_max))
        test_mod = SystModel(self, series, vars, z0=0.0, add=False)
        self._spec._shift_rf(0.0)
        xm = test_mod._xf
        hlenm = len(xm)//2
        ym = test_mod.eval(x=xm, params=test_mod._pars)
        ym_0 = np.ones(len(xm))
        ym_1 = np.concatenate([ym[:-hlenm], np.ones(hlenm)])
        ym_2 = np.concatenate([np.ones(hlenm), ym[hlenm:]])
        z_true = []
        chi2a = []
        chi2a_0 = []
        for i, z in enumerate(z_range):
            if verb:
                print(prefix, "I'm testing a %s system (logN=%2.2f, b=%2.2f) "
                      "at redshift %2.4f (%i/%i)…" \
                      % (series, logN, b, z, i+1, len(z_range)), end='\r')
            spec_temp._shift_rf(z)
            cond, chi2, chi2_0 = self._test(spec, xm, ym, ym_0, ym_1, ym_2, col,
                                            chi2_fact)
            if cond:
                z_true.append(z)
            chi2a.append(chi2)
            chi2a_0.append(chi2_0)
        spec_temp._shift_rf(0)

        #plt.plot(z_range, np.log(chi2a))
        #plt.plot(z_range, np.log(chi2a_0))
        #plt.show()
        if verb:
            print(prefix, "I've tested a %s system (logN=%2.2f, b=%2.2f) "\
                  "between redshifts %2.4f and %2.4f and found %i coincidences."\
                  % (series, logN, b, z_range[0], z_range[-1], len(z_true)))
        #self._update_spec()
        if len(z_true) > 0:
            self._add_fit(series, np.array(z_true), logN, b, resol,
                          chi2r_thres=chi2r_thres, verb=True)
        return 0

    def _update_spec(self):
        """ @brief Update spectrum with model """

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
        for i, r in enumerate(self._mods._t):
            mod = r['mod']
            model[s] = mod.eval(x=self._xs, params=mod._pars) * model[s]
        deabs[s] = cont[s] + y[s] - model[s]

        plt.plot(spec.x, spec._t['model'])
        plt.show()
