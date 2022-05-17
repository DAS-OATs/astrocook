from .functions import to_x, trans_parse
from .syst_model import SystModel
from astropy import table as at
from astropy import units as au
from copy import deepcopy as dc
import inspect
import logging
import numpy as np
import os
import pickle
from scipy.signal import find_peaks


class Feat():
    """Class for features"""


    def __init__(self, chunk):
        self._chunk = chunk
        self._left = chunk['x', 'model', 'cont'][0]
        self._right = chunk['x', 'model', 'cont'][-1]
        self._xunit = chunk['x'].unit
        self._yunit = chunk['y'].unit
        self._argmin = np.argmin(chunk['model']/chunk['cont'])


    def _ccf_compute(self, xmean, deltav):
        self._ccf_xmean = xmean
        self._ccf_deltav = deltav


    def _ew_compute(self):
        t = self._chunk
        self._ew = np.sum((t['xmax']-t['xmin']).to(au.nm)*(1-np.array(t['model']/t['cont'])))


    def _fwhm_compute(self):
        t = self._chunk
        hm = 1-0.5*(1-t['model'][self._argmin]/t['cont'][self._argmin])
        sel = np.where(t['model']/t['cont']<hm)
        self._fwhm = (t['x'].to(au.nm)[sel][-1]-t['x'].to(au.nm)[sel][0])


    def _logN_compute(self):
        if not self._systs_check(): return 0
        logN = np.array([s._pars['logN'] for s in self._systs.values()])
        dlogN = np.array([s._pars['dlogN'] for s in self._systs.values()])
        N = 10**logN
        dN = N*dlogN
        self._logN = np.log10(np.sum(N))
        self._dlogN = np.sqrt(np.sum(dN**2))/np.sum(N)


    def _snr_compute(self):
        t = self._chunk
        self._snr = np.nanmean(t['y']/t['dy'])


    def _systs_check(self):
        if not hasattr(self, '_systs'):
            logging.error("I cannot run %s before _systs_join."\
                          % inspect.stack()[1][3])
            return 0
        else:
            return 1


    def _systs_join(self, systs):
        self._systs_orig = systs
        self._systs = {}
        self._trans = {}
        for i, s in systs._d.items():
            for t, x in s._x.items():
                v = x.to(self._xunit).value
                if v>self._left[0] and v<self._right[0]:
                    self._systs[i] = s
                    self._trans[i] = t


    def _systs_logN_tot(self):
        self._dlogN_tot = None
        if hasattr(self, '_N_tot_par'):
            value = self._N_tot_par.value
            stderr = self._N_tot_par.stderr
            self._logN_tot = np.log10(value)
            if self._N_tot_par.stderr is not None:
                self._dlogN_tot = 0.5*np.log10((value+stderr)/(value-stderr))
            #print(self._logN_tot, self._dlogN_tot, 'est')
        else:
            self._logN_tot = self._logN_tot_par.value
            self._dlogN_tot = self._logN_tot_par.stderr
            #print(self._logN_tot, self._dlogN_tot)


    def _systs_stats(self):
        self._logN_compute()
        self._xz_compute()
        self._ew_compute()
        self._fwhm_compute()
        self._snr_compute()

    def _xz_compute(self):
        if not self._systs_check(): return 0
        w = [s._pars['logN']/s._pars['dz']**2 \
             if s._pars['dz']!=0 and not np.isnan(s._pars['dz']) \
             else s._pars['logN'] for s in self._systs.values()]
        z = [s._pars['z'] for s in self._systs.values()]
        #print(z)
        x = [to_x(zi, t).value for (zi, t) in zip(z, self._trans.values())]
        if np.nansum(w)==0:
            logging.warning("I cannot weight the positions of the systems "
                            "around z = %.3f by their errors and column "
                            "densities." % np.mean(z))
            self._z = np.mean(z)
            self._x = np.mean(x)*au.nm
        else:
            self._z = np.average(z, weights=w)
            self._x = np.average(x, weights=w)*au.nm



class FeatList(object):
    """Class for feature lists

    A FeatureList is a Frame with methods for handling absorption features."""

    def __init__(self):
        self._l = []
        #print('l', len(self._l))
        self._table_update()


    @property
    def t(self):
        return self._t


    def _add(self, chunk, systs):
        self._xunit = chunk['x'].unit
        self._yunit = chunk['y'].unit
        feat = Feat(chunk)
        feat._systs_join(systs)
        feat._systs_stats()
        self._l.append(feat)
        return 0


    def _check_attr(self, attr):
        if hasattr(self, attr):
            return getattr(self, attr)
        else:
            return None


    def create(self, spec, systs, thres, height=1e-1, prominence=1e-1):
        self._maxs_from_spec(spec, height, prominence)

        for i, f in enumerate(self._maxs[:-1]):
            fe = self._maxs[i+1]
            sel = np.s_[f:fe]
            modeln = spec._t['model'][sel]/spec._t['cont'][sel]
            min = np.amin(modeln)
            cut = np.where(modeln<1-thres*(1-min))
            if len(cut[0])>0:
                #plt.plot(spec._t[sel][cut]['x'], spec._t[sel][cut]['model'])
                self._add(spec._t[sel][cut], systs)

        self._table_update()
        #plt.show()




    def _load(self, new_dir, **kwargs):
        for file in sorted(os.listdir(new_dir)):
            with open(new_dir+file, 'rb') as f:
                o = pickle.load(f)
                if 'systs' in kwargs:
                    o._systs_join(kwargs['systs'])
                self._l.append(o)
        self._table_update()


    def _logN_tot(self, systs, set_specs=True, sel=None):
        systs._N_tot_specs = {}
        done = []
        if sel==None:
            l = self._l
        else:
            if type(sel)==str:
                sel_t = []
                for subs in sel.split(','):
                    subs_s = subs.split(':')
                    start = min(len(self._l), int(subs_s[0]))
                    end = min(len(self._l), int(subs_s[1]))
                    sel_t += range(start, end)
                sel = sel_t
            l = [self._l[s] for s in sel]
        print(len(l))
        for i, f in enumerate(l):
            s = max(f._systs.keys())
            syst = f._systs[s]
            if set_specs:
                if len(f._systs.keys())>1:
                    if s not in done:
                        lines_pref = 'lines_voigt_%s_' % str(s)
                        pars = dc(syst._mod._pars)
                        pdel = []
                        for p in pars:
                            if int(p.split('_')[2]) not in f._systs.keys(): pdel.append(p)
                        for p in pdel:
                            del pars[p]
                        N_tot = np.sum([10**pars[p] for p in pars if 'logN' in p])
                        N_other_expr = ''
                        for p in pars:
                            if 'logN' in p and 'logN_' not in p and p != lines_pref+'logN':
                                N_other_expr += '+10**%s' % p
                        systs._N_tot_specs[s] = (lines_pref, N_tot, N_other_expr)
                        done.append(s)
                        #pars.pretty_print()
            else:
                try:
                    f._N_tot_par = syst._mod._pars['lines_voigt_%s_N_tot' % str(s)]
                except:
                    f._logN_tot_par = syst._mod._pars['lines_voigt_%s_logN' % str(s)]
                f._systs_logN_tot()


    def _maxs_from_spec(self, spec, height=1e-1, prominence=1e-1):
        s = spec._where_safe
        maxs = find_peaks(spec._t['model'][s]/spec._t['cont'][s],
                          height=height, prominence=prominence)[0]
        self._maxs = np.hstack(([0], maxs, [-1]))


    def _save(self, new_dir):
        for i, o in enumerate(self._l):
            o._systs_orig = None
            o._systs = None
            with open(new_dir+'%04i.dat' % i, 'wb') as f:
                pickle.dump(o, f, pickle.HIGHEST_PROTOCOL)


    def _select_isolated(self, thres=1e-1):
        norm = self._t['y']/self._t['cont']
        sel = np.where(norm>1-thres)


    def _systs_update(self, systs):
        for i, f in enumerate(self._l):
            f._systs_join(systs)
            f._systs_stats()


    def _table_update(self):
        xl, x, m, c, r = [], [], [], [], []
        for f in self._l:
            xl.append(f._left[0])
            x.append(f._left[0])
            x.append(f._right[0])
            m.append(f._left[1])
            m.append(f._right[1])
            c.append(f._left[2])
            c.append(f._right[2])
        self._xleft = xl
        #self._t = FeatTable(x, m, c, self._check_attr('_xunit'),
        #                    self._check_attr('_yunit'))
        self._t = at.Table()
        self._t['x'] = at.Column(np.array(x, ndmin=1), dtype=float, unit=self._check_attr('_xunit'))
        self._t['model'] = at.Column(np.array(m, ndmin=1), dtype=float, unit=self._check_attr('_yunit'))
        self._t['cont'] = at.Column(np.array(c, ndmin=1), dtype=float, unit=self._check_attr('_yunit'))
        #print(self)
        #print(self._t)
        #for f in self._l:
        #    print(f._left[0], f._right[0])
        #print(len(self._t))

    def _z_lock(self):
        for i, f in enumerate(self._l):
            first = True
            for s in f._systs:
                #print(s)
                syst = f._systs[s]
                #print(f._systs[s]._pars)
                pars = syst._mod._pars
                z = 'lines_'+syst._func+'_'+str(s)+'_z'
                b = 'lines_'+syst._func+'_'+str(s)+'_b'
                #print(z)
                if first:
                    zref = z
                    first = False
                    dict = {b: (s, 'vary', False)}
                elif zref in pars:
                    #pars.pretty_print()
                    #pars[z].set(expr = zref+'+%.14f' % (pars[z]-pars[zref]))
                    dict = {z: (s, 'expr', zref+'+%.14f' % (pars[z]-pars[zref])),
                            b: (s, 'vary', False)}
                f._systs_orig._constrain(dict)
                    #print(f._systs_orig._constr)
                    #pars.pretty_print()
                #print(z, zref)
        #print(f._systs_orig._constr)

class FeatTable(at.Table):

    def __init__(self,
                 x=[],
                 model=[],
                 cont=[],
                 xunit=None,
                 yunit=None,
                 dtype=float):

        #super(FeatTable, self).__init__()
        self['x'] = at.Column(np.array(x, ndmin=1), dtype=dtype, unit=xunit)
        self['model'] = at.Column(np.array(model, ndmin=1), dtype=dtype, unit=yunit)
        self['cont'] = at.Column(np.array(cont, ndmin=1), dtype=dtype, unit=yunit)
