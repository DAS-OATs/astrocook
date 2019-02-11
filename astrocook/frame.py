from astropy import units as au
from astropy import table as at
from copy import deepcopy as dc
#from astropy.units import au.Quantity
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
        self._xunit = xunit
        self._yunit = yunit
        self._meta = meta
        self._dtype = dtype

    @property
    def t(self):
        return self._t

    @property
    def x(self):
        return au.Quantity(self._t['x'])

    @property
    def xmin(self):
        return au.Quantity(self._t['xmin'])

    @property
    def xmax(self):
        return au.Quantity(self._t['xmax'])

    @property
    def y(self):
        return au.Quantity(self._t['y'])

    @property
    def dy(self):
        return au.Quantity(self._t['dy'])

    @x.setter
    def x(self, val, dtype=float):
        self._t['x'] = np.array(val, dtype=dtype)
        self._t['x'].unit = val.unit

    @xmin.setter
    def xmin(self, val, dtype=float):
        self._t['xmin'] = np.array(val, dtype=dtype)
        self._t['xmin'].unit = val.unit

    @xmax.setter
    def xmax(self, val, dtype=float):
        self._t['xmax'] = np.array(val, dtype=dtype)
        self._t['xmax'].unit = val.unit

    @y.setter
    def y(self, val, dtype=float):
        self._t['y'] = np.array(val, dtype=dtype)
        self._t['y'].unit = val.unit

    @dy.setter
    def dy(self, val, dtype=float):
        self._t['dy'] = np.array(val, dtype=dtype)
        self._t['dy'].unit = val.unit

    @property
    def meta(self):
        return self._meta

    @meta.setter
    def meta(self, key, val):
        self._meta[key] = val

    def _safe(self, col):
        if isinstance(col, at.Column):
            col = au.Quantity(col)
        self._where_safe = ~np.isnan(col.value)
        return col[self._where_safe]

    def _copy(self, sel=None):
        if sel is None:
            sel = range(len(self.t))
        x = dc(self.x[sel])
        xmin = dc(self.xmin[sel])
        xmax = dc(self.xmax[sel])
        y = dc(self.y[sel])
        dy = dc(self.dy[sel])
        xunit = self._xunit
        yunit = self._yunit
        meta = self._meta
        dtype = self._dtype
        return type(self)(x, xmin, xmax, y, dy, xunit, yunit, meta, dtype)

    def _extract_region(self, xmin, xmax):
        """ @brief Extract a spectral region as a new frame.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return Spectral region
        """

        reg = dc(self)
        reg.x.unit.to(au.nm)
        where = np.full(len(reg.x), True)
        s = np.where(np.logical_and(self._safe(reg.x) > xmin,
                                    self._safe(reg.x) < xmax))
        where[s] = False
        reg._t.remove_rows(where)

        if len(reg.t) == 0:
            print(prefix, msg_output_fail)
            return None
        else:
            return reg
