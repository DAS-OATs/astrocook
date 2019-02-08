from astropy import units as au
from astropy import table as at
import numpy as np

class Frame(object):
    """Class for frames.

    A Frame is an astropy Table with the following columns:
        -# @x: channels;
        -# @xmin: lower limit for each channel;
        -# @xmax: upper limit for each channel;
        -# @y: flux density in the channel;
        -# @dy: error on @y.
    """

    def __init__(self,
                 x=[],
                 xmin=[],
                 xmax=[],
                 y=[],
                 dy=[],
                 xunit=au.nm,
                 yunit=au.erg/au.cm**2/au.s/au.nm,
                 meta={},
                 dtype=float):

        t = at.Table()
        t['x']  = at.Column(np.array(x, ndmin=1), dtype=dtype, unit=xunit)
        t['xmin'] = at.Column(np.array(xmin, ndmin=1), dtype=dtype, unit=xunit)
        t['xmax'] = at.Column(np.array(xmax, ndmin=1), dtype=dtype, unit=xunit)
        t['y']  = at.Column(np.array(y, ndmin=1) , dtype=dtype, unit=yunit)
        t['dy'] = at.Column(np.array(dy, ndmin=1), dtype=dtype, unit=yunit)
        self._t = t
        self._meta=meta

    @property
    def t(self):
        return self._t

    @property
    def x(self):
        return self._t['x']

    @property
    def xmin(self):
        return self._t['xmin']

    @property
    def xmax(self):
        return self._t['xmax']

    @property
    def y(self):
        return self._t['y']

    @property
    def dy(self):
        return self._t['dy']

    @x.setter
    def x(self, val, dtype=float):
        self._t['x'] = np.array(val, dtype=dtype)

    @xmin.setter
    def xmin(self, val, dtype=float):
        self._t['xmin'] = np.array(val, dtype=dtype)

    @xmax.setter
    def xmax(self, val, dtype=float):
        self._t['xmax'] = np.array(val, dtype=dtype)

    @y.setter
    def y(self, val, dtype=float):
        self._t['y'] = np.array(val, dtype=dtype)

    @dy.setter
    def dy(self, val, dtype=float):
        self._t['dy'] = np.array(val, dtype=dtype)

    @property
    def meta(self):
        return self._meta

    @meta.setter
    def meta(self, key, val):
        self._meta[key] = val
