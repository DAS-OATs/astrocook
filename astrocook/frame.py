from .message import *
from astropy import units as au
from astropy import constants as aconst
from astropy import table as at
from copy import deepcopy as dc
#from astropy.units import au.Quantity
import logging
import numpy as np

class Frame():
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
        self._rfz = 0.0

        self.x = au.Quantity(self._t['x'])

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

    @property
    def meta(self):
        return self._meta

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

    @meta.setter
    def meta(self, key, val):
        self._meta[key] = val


    def _append(self, frame):
        vstack = at.vstack([self._t, frame._t])
        if len(self._t) > 0:
            self._t = at.unique(vstack, keys=['x'])
        return 0


    def _copy(self, sel=None):
        """ @brief Copy a selection from a frame into a new frame.
        @param sel Selected rows. If 'None', all frame is copied.
        @return Copied frame
        """
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


    def _region_extract(self, xmin, xmax):
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
            logging.error(msg_output_fail)
            return None
        else:
            return reg


    def _safe(self, col):
        if isinstance(col, at.Column):
            col = au.Quantity(col)
        self._where_safe = ~np.isnan(col.value)
        return col[self._where_safe]


    def _shift_rf(self, z):
        """ @brief Shift to and from rest frame.
        @param z Redshift to use for shifting
        @return 0
        """

        fact = (1+self._rfz)/(1+z)
        self.x = self.x*fact
        self.xmin = self.xmin*fact
        self.xmax = self.xmax*fact
        self._rfz = z
        return 0


    def _shift_bary(self, v):
        """ @brief Shift to and from barycentic frame.
        @param v Velocity in the barycentric frame (km/s)
        @return 0
        """

        fact = 1+v/aconst.c.to(au.km/au.s).value
        self.x = self.x*fact
        self.xmin = self.xmin*fact
        self.xmax = self.xmax*fact
        return 0

    def _x_convert(self, zem=0, xunit=au.km/au.s):

        self._zem = zem
        xem = (1+zem) * 121.567*au.nm
        equiv = [(au.nm, au.km/au.s,
                  lambda x: np.log(x/xem.value)*aconst.c.to(au.km/au.s),
                  lambda x: np.exp(x/aconst.c.to(au.km/au.s).value)*xem.value)]

        self._xunit = xunit
        self._xunit_old = self.x.unit
        self.x = self.x.to(xunit, equivalencies=equiv)
        self.xmin = self.xmin.to(xunit, equivalencies=equiv)
        self.xmax = self.xmax.to(xunit, equivalencies=equiv)
        return 0


    def _y_convert(self, e_to_flux=None, yunit=au.erg/au.cm**2/au.s/au.nm):
        """ @brief Convert the y axis to electron or flux density units.
        @param e_to_flux Flux calibration array (same length as the spectrum)
        @param yunit Unit of electron or flux density
        @return 0
        """

        if e_to_flux == None:
            e_to_flux = np.ones(len(self.t))
        equiv = [(au.erg/au.cm**2/au.s/au.nm, au.electron/au.nm,
                  lambda y: y*e_to_flux, lambda y: y/e_to_flux)]

        self._yunit = yunit
        self.y = self.y.to(yunit, equivalencies=equiv)
        self.dy = self.dy.to(yunit, equivalencies=equiv)
        return 0

    def _y_scale(self, fact):
        self.y = self.y * fact
        self.dy = self.dy * fact
        return 0
