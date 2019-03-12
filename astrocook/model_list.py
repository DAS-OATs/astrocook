from astropy import table as at
import numpy as np

class ModelList(object):
    """ Class for model lists

    A SystList is a list of models. """

    def __init__(self,
                 sess,
                 mod=[],
                 chi2r=[],
#                 systs=[],
                 dtype=float):

        t = at.Table()
        t['mod'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        t['chi2r'] = at.Column(np.array(chi2r, ndmin=1), dtype=dtype)
#        t['systs'] = at.Column(np.array(mod, ndmin=1), dtype=object)
        self._t = t
        self._dtype = dtype

    @property
    def t(self):
        return self._t
