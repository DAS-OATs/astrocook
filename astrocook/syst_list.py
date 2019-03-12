#from .model import Model, ModelLines, ModelPSF
from .syst_model import SystModel
from .spectrum import Spectrum
from .vars import *
from .functions import convolve, lines_voigt, psf_gauss, running_mean
from astropy import table as at
from astropy import units as au
#from astropy import constants as ac
#from astropy.stats import sigma_clip
#from collections import OrderedDict
#from copy import deepcopy as dc
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
        t['func'] = at.Column(np.array(func, ndmin=1), dtype='S5')
        t['series'] = at.Column(np.array(func, ndmin=1), dtype='S100')
        t['z'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        self._t = t
        self._t['mod'] = np.empty(len(self.z), dtype=object)
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
        for i, r in enumerate(self._t):
            mod = r['mod']
            model[s] = mod.eval(x=self._xs, params=mod._pars) * model[s]
        deabs[s] = cont[s] + y[s] - model[s]

    def add(self, series='Ly_a', z=2.0, logN=13, b=10, resol=35000):
        """ @brief Add a Voigt model for a system.
        @param series Series of transitions
        @param z Redshift
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

        mod = SystModel(self._spec, series, vars)
        self._t.add_row(['voigt_func', series, z, mod, 0.0])
        self._update_spec()

        return 0

    def fit(self, row=0):
        """ @brief Fit a Voigt model for a system.
        @param row Row in the system list (0 means last row)
        @return 0
        """

        row = int(row)

        mod = self._t[row]['mod']
        mod = mod.fit()

        return 0
