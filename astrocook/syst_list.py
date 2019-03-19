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
                 #spec=None,
                 sess=None,
                 series=[],
                 func=[],
                 z=[],
                 logN=[],
                 b=[],
                 mod=[],
                 chi2r=[],
                 id=[],
                 dtype=float):

        self._spec = sess.spec
        #self._mods = ModelList()#sess)
        #self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))  # Full x array, without NaNs
        #self._s = self._spec._where_safe
        self._id = 0

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
        self._t['id'] = np.empty(len(self.z), dtype=int)

        self._mods_t = at.Table()
        self._mods_t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype)
        self._mods_t['mod'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        self._mods_t['chi2r'] = at.Column(np.array(chi2r, ndmin=1), dtype=dtype)
        self._mods_t['id'] = at.Column(np.array(id, ndmin=1), dtype=object)

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

    """
    def _add(self, series='Ly_a', z=2.0, logN=13, b=10, resol=70000):

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        vars = {'z': z, 'logN': logN, 'b': b, 'resol': resol}
        id = len(self._t)

        self._mod = SystModel(self, series, vars, z0=z)
        self._t.add_row(['voigt_func', series, z, z, logN, b, None, self._id])

        return 0
    """

    def _add2(self, series='Ly_a', z=2.0, logN=13, b=10, resol=70000):
        """ @brief Add a system to a system list.
        """

        self._t.add_row(['voigt_func', series, z, z, logN, b, None, self._id])

        return 0

    """
    def _add_fit(self, series='Ly_a', z=2.0, logN=14, b=10, resol=70000,
                 chi2r_thres=None, fit_kws={}, verb=False):

        z_range = z if np.size(z) > 1 else [z]
        for i, z in enumerate(z_range):
            if verb:
                print(prefix, "I'm fitting a %s system at redshift %2.4f "\
                      "(%i/%i)…" % (series, z, i+1, len(z_range)), end='\r')
            self._add(series, z, logN, b, resol)
            self._fit(chi2r_thres, fit_kws)
        if chi2r_thres != np.inf:
            z_rem = self._clean(chi2r_thres)
        else:
            z_rem = []

        if verb:
            print(prefix, "I've fitted %i %s systems between redshift %2.4f "\
                  "and %2.4f." % (len(z_range), series, z_range[0], z_range[-1]))
        return 0
    """

    def _append(self, frame):
        vstack_t = at.vstack([self._t, frame._t])
        vstack_mods_t = at.vstack([self._mods_t, frame._mods_t])
        self._t = at.unique(vstack_t, keys=['z0'])
        self._mods_t = at.unique(vstack_mods_t, keys=['z0'])
        return 0

    def _clean(self, chi2r_thres=np.inf):

        rem = np.where(self._t['chi2r'] > chi2r_thres)[0]
        #mods_rem = np.where(self._mods._t['chi2r'] > chi2r_thres)[0]
        mods_rem = np.where(self._mods_t['chi2r'] > chi2r_thres)[0]
        z_rem = self._t['z'][rem]
        if rem != []:
            self._t.remove_rows(rem)
        if mods_rem != []:
            self._mods_t.remove_rows(mods_rem)
        return 0

    """
    def _fit(self, chi2r_thres=None, fit_kws={}):

        if chi2r_thres == None:
            chi2r_thres = np.inf

        self._mod._fit(fit_kws)

        #mod = np.where(self._mod == self._mods._t['mod'])[0][0]
        mod = np.where(self._mod == self._mods_t['mod'])[0][0]
        #id = self._mods._t['id'][mod]
        id = self._mods_t['id'][mod]
        for i in id:
            iw = np.where(self._t['id']==i)[0][0]
            pref = 'lines_voigt_'+str(i)
            self._t[iw]['z'] = self._mod._pars[pref+'_z'].value
            self._t[iw]['logN'] = self._mod._pars[pref+'_logN'].value
            self._t[iw]['b'] = self._mod._pars[pref+'_b'].value
            self._t[iw]['chi2r'] = self._mod._chi2r
        self._id += 1

        return 0
    """

    def _update(self, mod):

        if mod._group_sel == -1:
            self._mods_t.add_row([mod._z0, mod, None, []])
        else:
            self._mods_t[mod._group_sel]['mod'] = mod
        self._mods_t[mod._group_sel]['chi2r'] = mod._chi2r
        self._mods_t[mod._group_sel]['id'].append(mod._id)


        modw = np.where(mod == self._mods_t['mod'])[0][0]
        id = self._mods_t['id'][modw]
        for i in id:
            iw = np.where(self._t['id']==i)[0][0]
            pref = 'lines_voigt_'+str(i)
            self._t[iw]['z'] = mod._pars[pref+'_z'].value
            self._t[iw]['logN'] = mod._pars[pref+'_logN'].value
            self._t[iw]['b'] = mod._pars[pref+'_b'].value
            self._t[iw]['chi2r'] = mod._chi2r
        self._id += 1

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

    """
    def _test_fit(self, spec, series='Ly_a', z=2.0, logN=14, b=10, resol=70000,
                  col='y', chi2_fact=1.0, chi2r_thres=2.0, fit_kws={},
                  verb=False):

        z_range = z if np.size(z) > 1 else [z]

        #spec_temp = spec
        vars = {'z': 0.0, 'logN': logN, 'b': b, 'resol': resol}

        z_min = np.min([(np.min(spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_max = np.max([(np.max(spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])

        z_range = z_range[np.where(np.logical_and(z_range > z_min,
                                                  z_range < z_max))]
        spec._shift_rf(0.5*(z_min+z_max))
        test_mod = SystModel(self, series, vars, z0=0.0, add=False)
        spec._shift_rf(0.0)
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
            spec._shift_rf(z)
            cond, chi2, chi2_0 = self._test(spec, xm, ym, ym_0, ym_1, ym_2, col,
                                            chi2_fact)
            if cond:
                z_true.append(z)
            chi2a.append(chi2)
            chi2a_0.append(chi2_0)
        spec._shift_rf(0)

        if verb:
            print(prefix, "I've tested a %s system (logN=%2.2f, b=%2.2f) "\
                  "between redshifts %2.4f and %2.4f and found %i coincidences."\
                  % (series, logN, b, z_range[0], z_range[-1], len(z_true)))

        if len(z_true) > 0:
            self._add_fit(series, np.array(z_true), logN, b, resol,
                          chi2r_thres=chi2r_thres, fit_kws=fit_kws,
                          verb=True)

        return 0
    """


class SystList2(SystList):
    """ Class for system lists

    A SystList is a list of absorption systems with methods for handling
    spectral lines. """

    def __init__(self,
                 series=[],
                 func=[],
                 z=[],
                 logN=[],
                 b=[],
                 mod=[],
                 chi2r=[],
                 id=[],
                 dtype=float):

        self._id = 0

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
        self._t['chi2r'] = np.empty(len(self.z), dtype=dtype)
        self._t['id'] = np.empty(len(self.z), dtype=int)

        self._mods_t = at.Table()
        self._mods_t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype)
        self._mods_t['mod'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        self._mods_t['chi2r'] = at.Column(np.array(chi2r, ndmin=1), dtype=dtype)
        self._mods_t['id'] = at.Column(np.array(id, ndmin=1), dtype=object)

        self._dtype = dtype
