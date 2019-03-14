from astropy import table as at
import numpy as np

class ModelList(object):
    """ Class for model lists

    A SystList is a list of models. """

    def __init__(self,
                 #sess,
                 z=[],
                 mod=[],
                 chi2r=[],
#                 systs=[],
                 dtype=float):

        t = at.Table()
        t['z0'] = at.Column(np.array(z, ndmin=1), dtype=dtype)
        t['mod'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        t['chi2r'] = at.Column(np.array(chi2r, ndmin=1), dtype=dtype)
#        t['systs'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        self._t = t
        self._dtype = dtype

    @property
    def t(self):
        return self._t

    def _append(self, frame):
        vstack = at.vstack([self._t, frame._t])
        self._t = at.unique(vstack, keys=['z0'])
        return 0
