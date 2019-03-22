from .vars import *
from .functions import convolve, lines_voigt, psf_gauss, running_mean
from astropy import table as at
from astropy import units as au
#from matplotlib import pyplot as plt
from copy import deepcopy as dc
import numpy as np

prefix = "System list:"

class SystList(object):
    """ Class for system lists

    A SystList is a list of absorption systems with methods for handling
    spectral lines. """

    def __init__(self,
                 id_start=0,
                 series=[],
                 func=[],
                 z=[],
                 logN=[],
                 b=[],
                 mod=[],
                 chi2r=[],
                 id=[],
                 dtype=float):

        self._id = id_start

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

        mods_t = at.Table()
        mods_t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype)
        mods_t['mod'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        #mods_t['chi2r'] = at.Column(np.array(chi2r, ndmin=1), dtype=dtype)
        #mods_t['id'] = at.Column(np.array(id, ndmin=1), dtype=object)
        self._mods_t = mods_t
        self._mods_t['chi2r'] = np.empty(len(self.z), dtype=dtype)
        self._mods_t['id'] = np.empty(len(self.z), dtype=object)

        self._dtype = dtype


    @property
    def t(self):
        return self._t

    @property
    def mods_t(self):
        return self._mods_t

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


    def _append(self, frame, unique=True):
        vstack_t = at.vstack([self._t, frame._t])
        vstack_mods_t = at.vstack([self._mods_t, frame._mods_t])
        if unique:
            self._t = at.unique(vstack_t, keys=['z0', 'z'])
            self._mods_t = at.unique(vstack_mods_t, keys=['z0'])
        return 0


    def _clean(self, chi2r_thres=np.inf, verb=True):

        rem = np.where(self._t['chi2r'] > chi2r_thres)[0]
        #mods_rem = np.where(self._mods._t['chi2r'] > chi2r_thres)[0]
        mods_rem = np.where(self._mods_t['chi2r'] > chi2r_thres)[0]
        #print(rem, mods_rem)
        z_rem = self._t['z'][rem]
        if len(rem) > 0:
            self._t.remove_rows(rem)
            print(prefix, "I removed %i systems because they had a "\
                  "chi-squared above %2.2f." % (len(rem), chi2r_thres))
        #print(rem, mods_rem)
        if len(mods_rem) > 0:
            #print("ecco")
            self._mods_t.remove_rows(mods_rem)
        #print(self._mods_t['z0', 'chi2r', 'id'])
        return 0

    def _freeze(self):
        """ Create a frozen copy of the tables self._t and self._mods_t
        """
        t = dc(self._t)
        mods_t = dc(self._mods_t)

        # Needed, otherwise the id objects are not copied
        for i, m in enumerate(mods_t):
            m['id'] = dc(self._mods_t['id'][i])

        return t, mods_t

    def _unfreeze(self, t, mods_t):
        """ Restore from a frozen copy of the tables self._t and self._mods_t
        """

        self._t = t
        self._mods_t = mods_t


    def _update(self, mod):

        #print("inside")
        if mod._group_sel == -1:
            self._mods_t.add_row([mod._z0, mod, None, []])
        else:
            self._mods_t[mod._group_sel]['mod'] = mod
        self._mods_t[mod._group_sel]['chi2r'] = mod._chi2r
        """
        try:
            print(sess._mods_t_old['chi2r','id'])
            print(hex(id(sess._mods_t_old[mod._group_sel]['id'])))
            print(hex(id(self._mods_t[mod._group_sel]['id'])))
            print(hex(id(sess._mods_t_old['id'])))
            print(hex(id(self._mods_t['id'])))
        except:
            pass
        """
        self._mods_t[mod._group_sel]['id'].append(mod._id)

        """
        try:
            print(sess._mods_t_old['chi2r','id'])
        except:
            pass
        """
        modw = np.where(mod == self._mods_t['mod'])[0][0]
        ids = self._mods_t['id'][modw]
        for i in ids:
            iw = np.where(self._t['id']==i)[0][0]
            pref = 'lines_voigt_'+str(i)
            self._t[iw]['z'] = mod._pars[pref+'_z'].value
            self._t[iw]['logN'] = mod._pars[pref+'_logN'].value
            self._t[iw]['b'] = mod._pars[pref+'_b'].value
            self._t[iw]['chi2r'] = mod._chi2r
        self._id += 1
        #print("until here")
