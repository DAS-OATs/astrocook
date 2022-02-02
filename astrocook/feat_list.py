from .functions import to_x, trans_parse
from astropy import table as at
from astropy import units as au
import numpy as np
import pickle
from scipy.signal import argrelmax


class Feat():
    """Class for features"""


    def __init__(self, chunk):
        self._chunk = chunk
        self._left = chunk['x', 'model', 'cont'][0]
        self._right = chunk['x', 'model', 'cont'][-1]
        self._xunit = chunk['x'].unit
        self._yunit = chunk['y'].unit


    def _systs_join(self, systs, feats):
        self._systs = {}
        self._trans = {}
        for i, s in systs._d.items():
            for t, x in s._x.items():
                v = x.to(self._xunit).value
                if v>self._left[0] and v<self._right[0]:
                    self._systs[i] = s
                    self._trans[i] = t


    def _systs_ave(self):
        N = [10**s._pars['logN'] for s in self._systs.values()]
        self._logN = np.mean(np.log10(np.sum(N)))
        w = [s._pars['logN']/s._pars['dz']**2 for s in self._systs.values()]
        z = [s._pars['z'] for s in self._systs.values()]
        x = [to_x(zi, t).value for (zi, t) in zip(z, self._trans.values())]
        self._z = np.average(z, weights=w)*au.nm
        self._x = np.average(x, weights=w)*au.nm
        #print(self._z, z)
        #print(self._x, x)


class FeatList(object):
    """Class for feature lists

    A FeatureList is a Frame with methods for handling absorption features."""

    def __init__(self,
                 l=[]):
        self._l = l
        self._table_update()


    @property
    def t(self):
        return self._t


    def _add(self, chunk, systs):
        self._xunit = chunk['x'].unit
        self._yunit = chunk['y'].unit
        feat = Feat(chunk)
        feat._systs_join(systs, self)
        feat._systs_ave()
        self._l.append(feat)
        return 0


    def _argrelmax_from_spec(self, spec):
        s = spec._where_safe
        armax = argrelmax(spec._t['model'][s]/spec._t['cont'][s])[0]
        self._argrelmax = np.hstack(([0], armax, [-1]))


    def _check_attr(self, attr):
        if hasattr(self, attr):
            return getattr(self, attr)
        else:
            return None


    def _open(self, new_dir):
        for file in os.listdir(new_dir):
            with open(new_dir+file, 'rb') as f:
                self._l.append(pickle.load(f))
        self._table_update()


    def _save(self, new_dir):
        for i, o in enumerate(self._l):
            with open(new_dir+'%04i.dat' % i, 'wb') as f:
                pickle.dump(o, f, pickle.HIGHEST_PROTOCOL)


    def _select_isolated(self, thres=1e-1):
        norm = self._t['y']/self._t['cont']
        sel = np.where(norm>1-thres)


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


    def create(self, spec, systs, thres):
        self._argrelmax_from_spec(spec)

        for i, f in enumerate(self._argrelmax[:-1]):
            fe = self._argrelmax[i+1]
            sel = np.s_[f:fe]
            cut = np.where(spec._t['model'][sel]/spec._t['cont'][sel]<1-thres)
            if len(cut[0])>0:
                self._add(spec._t[sel][cut], systs)

        self._table_update()


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
