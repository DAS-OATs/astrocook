from .vars import *
from .functions import convolve, lines_voigt, psf_gauss, running_mean
from astropy import table as at
from astropy import units as au
#from matplotlib import pyplot as plt
import numpy as np

prefix = "System list:"

class SystList(object):
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
        """ @brief Add a system to a system list.
        """

        self._t.add_row(['voigt_func', series, z, z, logN, b, None, self._id])

        return 0


    def _append(self, frame):
        vstack_t = at.vstack([self._t, frame._t])
        vstack_mods_t = at.vstack([self._mods_t, frame._mods_t])
        self._t = at.unique(vstack_t, keys=['z0'])
        self._mods_t = at.unique(vstack_mods_t, keys=['z0'])
        return 0


    def _clean(self, chi2r_thres=np.inf, verb=True):

        rem = np.where(self._t['chi2r'] > chi2r_thres)[0]
        #mods_rem = np.where(self._mods._t['chi2r'] > chi2r_thres)[0]
        mods_rem = np.where(self._mods_t['chi2r'] > chi2r_thres)[0]
        z_rem = self._t['z'][rem]
        if rem != []:
            self._t.remove_rows(rem)
            print(prefix, "I removed %i systems because they had a "\
                  "chi-squared below %2.2f." % (len(rem), chi2r_thres))
        if mods_rem != []:
            self._mods_t.remove_rows(mods_rem)
        return 0


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
