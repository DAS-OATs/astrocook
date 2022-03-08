from .vars import *
from .functions import convolve, lines_voigt, psf_gauss, running_mean, to_x, trans_parse
from .message import msg_output_fail
from astropy import table as at
from astropy import units as au
#from matplotlib import pyplot as plt
from copy import deepcopy as dc
import logging
import numpy as np

class SystList(object):
    """ Class for system lists

    A SystList is a list of absorption systems with methods for handling
    spectral lines. """

    def __init__(self,
                 id_start=0,
                 func=[],
                 series=[],
                 z=[],
                 dz=[],
                 logN=[],
                 dlogN=[],
                 b=[],
                 db=[],
                 btur=[],
                 dbtur=[],
                 mod=[],
                 resol=[],
                 chi2r=[],
                 snr=[],
                 id=[],
                 meta={},
                 dtype=float):

        self._id = id_start
        self._constr = {}

        t = at.Table()
        zunit = au.dimensionless_unscaled
        logNunit = au.dimensionless_unscaled
        bunit = au.km/au.s
        t['func'] = at.Column(np.array(func, ndmin=1), dtype='S5')
        t['series'] = at.Column(np.array(series, ndmin=1), dtype='S100')
        t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        t['z'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        t['dz'] = at.Column(np.array(dz, ndmin=1), dtype=dtype, unit=zunit)
        t['logN'] = at.Column(np.array(logN, ndmin=1), dtype=dtype,
                              unit=logNunit)
        t['dlogN'] = at.Column(np.array(dlogN, ndmin=1), dtype=dtype,
                               unit=logNunit)
        t['b'] = at.Column(np.array(b, ndmin=1), dtype=dtype, unit=bunit)
        t['db'] = at.Column(np.array(db, ndmin=1), dtype=dtype, unit=bunit)
        t['btur'] = at.Column(np.array(btur, ndmin=1), dtype=dtype, unit=bunit)
        t['dbtur'] = at.Column(np.array(dbtur, ndmin=1), dtype=dtype, unit=bunit)
        self._t = t
        #if resol != []:
        if len(resol)==len(self.z) and len(resol)>0:
            self._t['resol'] = resol
        else:
            self._t['resol'] = np.empty(len(self.z), dtype=dtype)
#        if chi2r != []:
        if len(chi2r)==len(self.z) and len(chi2r)>0:
            self._t['chi2r'] = chi2r
        else:
            self._t['chi2r'] = np.empty(len(self.z), dtype=dtype)
        if len(snr)==len(self.z) and len(snr)>0:
            self._t['snr'] = snr
        else:
            self._t['snr'] = np.empty(len(self.z), dtype=dtype)
            self._t['snr'] = np.nan
#        if id != []:
        if len(id)==len(self.z) and len(id)>0:
            self._t['id'] = id
        else:
            self._t['id'] = np.empty(len(self.z), dtype=int)
        mods_t = at.Table()
        mods_t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype)

        # Currently models cannot be saved, so they can't be retrieved from a
        # saved session. This 'try' is meant to skip model definition when a
        # session is loaded from a file.
        try:
            mods_t['mod'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        except:
            mods_t['mod'] = None
        #mods_t['chi2r'] = at.Column(np.array(chi2r, ndmin=1), dtype=dtype)
        #mods_t['id'] = at.Column(np.array(id, ndmin=1), dtype=object)
        self._mods_t = mods_t
        self._mods_t['chi2r'] = np.empty(len(self.z), dtype=dtype)
        self._mods_t['id'] = np.empty(len(self.z), dtype=object)

        self._meta = meta
        #self._meta = self._constr
        self._dtype = dtype

        self._compressed = False

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


    ### Deprecated, kept for backward compatibility with systs_new_from_resids
    def _add(self, series='Ly_a', z=2.0, logN=13, b=10, resol=70000):

        self._t.add_row(['voigt', series, z, z, None, logN, None, b, None,
                         None, None, self._id])

        return 0

    def _append(self, frame, unique=False):
        vstack_t = at.vstack([self._t, frame._t])
        vstack_mods_t = at.vstack([self._mods_t, frame._mods_t])

        #print(np.array(self._mods_t['z0']))
        #print(np.array(frame._mods_t['id']))
        #print(np.array(vstack_mods_t['id']))

        if unique:
            self._t = at.unique(vstack_t, keys=['z0', 'z'])
#            self._mods_t = at.unique(vstack_mods_t, keys=['z0'])
            self._mods_t = at.unique(vstack_mods_t, keys=['z0', 'id'])
        else:
            self._t = vstack_t
            self._mods_t = vstack_mods_t
        #print(self._mods_t['z0', 'id'])
        #print(len(self._mods_t))
        return 0


    def _clean(self, chi2r_thres=np.inf, verb=True):

        rem = np.where(self._t['chi2r'] > chi2r_thres)[0]
        mods_rem = np.where(self._mods_t['chi2r'] > chi2r_thres)[0]
        z_rem = self._t['z'][rem]
        if len(rem) > 0:
            self._t.remove_rows(rem)
            logging.info("I removed %i systems because they had a "\
                         "chi-squared above %2.2f." % (len(rem), chi2r_thres))
        if len(mods_rem) > 0:
            self._mods_t.remove_rows(mods_rem)
        return 0


    def _compress(self):
        if not self._compressed:
            self._t_uncompressed = dc(self._t)
            self._t['group'] = np.empty(len(self._t), dtype=int)
            for i, ids in enumerate(self._mods_t['id']):
                for id in ids:
                    self._t['group'][self._t['id']==id] = i
            t_by_group = self._t.group_by('group')
            t = at.Table()
            #t = t_by_group.groups.aggregate(np.mean)
            t['series'] = at.Column(np.array([g['series'][len(g)//2] \
                                              for g in t_by_group.groups]))
            t['z'] = at.Column(np.array([np.average(g['z'], weights=10**g['logN']) \
                                         for g in t_by_group.groups]))
            t['logN'] = at.Column(np.array([np.log10(np.sum(10**g['logN'])) \
                                            for g in t_by_group.groups]))
            for c in t_by_group.colnames[10:-1]:
                t[c] = at.Column(np.array([g[c][len(g)//2] \
                                           for g in t_by_group.groups]))
            """
            t['chi2r'] = at.Column(np.array([g['chi2r'][len(g)//2] \
                                             for g in t_by_group.groups]))
            t['id'] = at.Column(np.array([g['id'][len(g)//2] \
                                          for g in t_by_group.groups]))
            """
            t.sort(['z','id'])
            self._t = t
            self._compressed = True
        else:
            self._t = self._t_uncompressed
            self._compressed = False


    def _constrain(self, dict):
        #self._constr = {}
        #print(self)
        for k, v in dict.items():
            #print(k, dict[k])
            for m in self._mods_t:
                if v[0] in m['id']:
                    """
                    if k in ['lines_voigt_%i_b' %i for i in (33,34,35,40,41,42)]:
                        print('inside')
                        print(v)
                        m['mod']._pars.pretty_print()
                    """
                    if v[1]=='expr':
                        #print(k, v[2], type(v[2]))
                        m['mod']._pars[k].set(expr=v[2])
                        if v[2]=='':
                            m['mod']._pars[k].set(vary=True)
                        self._constr[k] = (v[0], k.split('_')[-1], v[2])
                    if v[1]=='vary':
                        m['mod']._pars[k].set(vary=v[2])
                        if v[2]:
                            if k in self._constr: del self._constr[k]
                        else:
                            self._constr[k] = (v[0], k.split('_')[-1], None)
                        #print(v[0], v[1], v[2])
                        #print(m['mod']._pars[k].__dict__)
                    #m['mod']._pars.pretty_print()
        #print(self._constr)
                #print(m['mod']._pars)
                #m['mod']._pars.pretty_print()
                #print('[',id(m),']')

    def _freeze(self):
        """ Create a frozen copy of the tables self._t and self._mods_t
        """
        t = dc(self._t)
        mods_t = dc(self._mods_t)

        # Needed, otherwise the id objects are not copied
        for i, m in enumerate(mods_t):
            m['id'] = dc(self._mods_t['id'][i])

        return t, mods_t


    def _region_extract(self, xmin, xmax):
        """ @brief Extract a spectral region as a new syst_list.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return Spectral region
        """

        reg = dc(self)
        #reg_x = np.array([to_x(z, series_d[s][0]).value for (z, s) in zip(reg._t['z0'], reg._t['series'])])*au.nm
        reg_x = np.array([to_x(z, trans_parse(s)[0]).value for (z, s) in zip(reg._t['z0'], reg._t['series'])])*au.nm
        where = np.full(len(reg_x), True)
        s = np.where(np.logical_and(reg_x > xmin, reg_x < xmax))
        where[s] = False
        reg._t.remove_rows(where)
        if len(reg.t) == 0:
            #logging.error(msg_output_fail)
            return None
        else:
            return reg


    def _unfreeze(self, t, mods_t):
        """ Restore from a frozen copy of the tables self._t and self._mods_t
        """

        self._t = t
        self._mods_t = mods_t


    def _update(self, mod, mod_t=True, t=True):

        #print(mod._id, mod._group_sel)
        if mod_t:
            if mod._group_sel == -1:
                self._mods_t.add_row([mod._z0, mod, None, []])
            else:
                self._mods_t[mod._group_sel]['mod'] = mod
            try:
                self._mods_t[mod._group_sel]['chi2r'] = mod._chi2r
            except:
                self._mods_t[mod._group_sel]['chi2r'] = np.nan
            #print(self._mods_t[mod._group_sel]['id'], mod._id)
            self._mods_t[mod._group_sel]['id'].append(mod._id)
            #print(self._mods_t[mod._group_sel]['id'])

        if t:
            modw = np.where(mod == self._mods_t['mod'])[0][0]
            ids = self._mods_t['id'][modw]
            #print(ids)
            for i in ids:
                try:
                    iw = np.where(self._t['id']==i)[0][0]
                    pref = 'lines_voigt_'+str(i)
                    self._t[iw]['z'] = mod._pars[pref+'_z'].value
                    self._t[iw]['dz'] = mod._pars[pref+'_z'].stderr
                    self._t[iw]['logN'] = mod._pars[pref+'_logN'].value
                    self._t[iw]['dlogN'] = mod._pars[pref+'_logN'].stderr
                    self._t[iw]['b'] = mod._pars[pref+'_b'].value
                    self._t[iw]['db'] = mod._pars[pref+'_b'].stderr
                    self._t[iw]['btur'] = mod._pars[pref+'_btur'].value
                    self._t[iw]['dbtur'] = mod._pars[pref+'_btur'].stderr
                    try:
                        self._t[iw]['chi2r'] = mod._chi2r
                    except:
                        self._t[iw]['chi2r'] = np.nan
                except:
                    pass


        self._id += 1

        #print(self._mods_t['id', 'chi2r'])
        #print(self._t)
