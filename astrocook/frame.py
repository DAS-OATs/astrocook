from .functions import x_convert
from .message import *
from .vars import inoue, logN_index, logN_lims, tau_index, tau_norm, xem_d
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
        self._meta = dc(meta)
        self._meta_backup = dc(meta)
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
        if self._xunit != frame._xunit:
            frame._t['x'] = frame._t['x'].to(self._xunit)
            frame._t['xmin'] = frame._t['xmin'].to(self._xunit)
            frame._t['xmax'] = frame._t['xmax'].to(self._xunit)
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

    def _lya_corr_basic(self, zem, logN_thres, input_col='y', apply=True,
                        verb=True):

        dv_prox = 1e4
        zprox = zem - (1.0+zem) * dv_prox / aconst.c.to(au.km/au.s).value

        num = 1
        den = 1
        fact_0 = 1
        fact_1 = 1e14 / np.log(1e14)
        fact_2 = np.log(1e18) / np.log(1e14) * 1e5;

        if logN_thres != 100:
            if logN_thres <= logN_lims[0]:
                num = 0
            if logN_thres > logN_lims[0] and logN_thres <= logN_lims[1]:
                num = (np.power(10, logN_thres*(2+logN_index)) \
                      - np.power(10, logN_lims[0] * (2+logN_index))) \
                      / (2.0+logN_index)
            if logN_thres > logN_lims[1] and logN_thres <= logN_lims[2]:
                num = (np.power(10, logN_lims[1] * (2+logN_index)) \
                      - np.power(10, logN_lims[0] * (2+logN_index))) \
                      / (2+logN_index)
                num = num + fact_1 \
                      * (np.power(10, logN_thres * (1+logN_index))) \
                        / (1+logN_index) \
                      * (np.log(np.power(10, logN_thres)) - 1/(1+logN_index))
                num = num - fact_1 \
                      * (np.power(10, logN_lims[1] * (1+logN_index))) \
                        / (1+logN_index) \
                      * (np.log(np.power(10, logN_lims[1])) - 1/(1+logN_index))
            if logN_thres > logN_lims[2] and logN_thres <= logN_lims[3]:
                num = (np.power(10, logN_lims[1] * (2+logN_index)) \
                      - np.power(10, logN_lims[0] * (2+logN_index))) \
                      / (2+logN_index)
                num = num + fact_1 \
                      * (np.power(10, logN_lims[2] * (1+logN_index))) \
                        / (1+logN_index) \
                      * (np.log(np.power(10, logN_lims[2])) - 1/(1+logN_index))
                num = num - fact_1 \
                      * (np.power(10, logN_lims[1] * (1+logN_index))) \
                        / (1+logN_index) \
                      * (np.log(np.power(10, logN_lims[1])) - 1/(1+logN_index))
                num = num + fact_2 \
                      * (np.power(10, logN_thres * (1.5 + logN_index)) \
                      - np.power(10, logN_lims[2] * (1.5 + logN_index))) \
                      / (1.5 + logN_index)

            if logN_thres > logN_lims[3]:
                num = 1;
            den = (np.power(10, logN_lims[1] * (2+logN_index)) \
                  - np.power(10, logN_lims[0] * (2+logN_index))) \
                  / (2+logN_index)
            den = den + fact_1 \
                  * (np.power(10, logN_lims[2] * (1+logN_index))) \
                    / (1+logN_index) \
                  * (np.log(np.power(10, logN_lims[2])) - 1/(1+logN_index))
            den = den - fact_1 \
                  * (np.power(10, logN_lims[1] * (1+logN_index))) \
                    / (1+logN_index)\
                  * (np.log(np.power(10, logN_lims[1])) - 1/(1+logN_index))

        corr = np.ones(len(self.y))
        z = (self.x.to(au.nm)/xem_d['Ly_a']).value - 1
        no_prox = z<zprox
        prox = np.logical_and(z>zprox,z<zem)
        corr[no_prox] = np.exp(tau_norm * np.power(1+z[no_prox], tau_index) \
                               * num/den) * fact_0
        corr[prox] = (1 - np.exp(tau_norm * np.power(1+zprox, tau_index) \
                                 * num / den) * fact_0) \
                     / (zem - zprox) * (z[prox] - zprox) \
                     + np.exp(tau_norm * np.power(1+zprox, tau_index) \
                              * num / den) * fact_0

        if apply:
            if verb:
                if input_col+'_taucorr' not in self._t.colnames:
                    logging.info("I'm adding column %s_taucorr." % input_col)
            taucorr = dc(self._t[input_col])
            taucorr *= corr
            self._t[input_col+'_taucorr'] = taucorr

        return corr

    def _lya_corr_inoue(self, zem, input_col='y', apply=True):

        x = self.x.to(au.nm).value

        x_ll = xem_d['Ly_lim'].to(au.nm).value

        tau_ls_laf = np.zeros(len(x))
        tau_ls_dla = np.zeros(len(x))
        tau_lc_laf = np.zeros(len(x))
        tau_lc_dla = np.zeros(len(x))

        for (i, xi, a_laf_1, a_laf_2, a_laf_3, a_dla_1, a_dla_2) in inoue:
            xi = 0.1*xi

            # Eq. 21
            s_1 = np.where(np.logical_and(x>xi, np.logical_and(x<xi*2.2, x<xi*(1+zem))))[0]
            s_2 = np.where(np.logical_and(x>xi, np.logical_and(np.logical_and(x>xi*2.2, x<xi*5.7), x<xi*(1+zem))))[0]
            s_3 = np.where(np.logical_and(x>xi, np.logical_and(np.logical_and(x>xi*2.2, x>xi*5.7), x<xi*(1+zem))))[0]
            tau_ls_laf[s_1] = tau_ls_laf[s_1] + a_laf_1 * (x[s_1]/xi)**1.2
            tau_ls_laf[s_2] = tau_ls_laf[s_2] + a_laf_2 * (x[s_2]/xi)**3.7
            tau_ls_laf[s_3] = tau_ls_laf[s_3] + a_laf_3 * (x[s_3]/xi)**5.5

            # Eq. 22
            s_4 = np.where(np.logical_and(x>xi, np.logical_and(x<xi*3.0, x<xi*(1+zem))))[0]
            s_5 = np.where(np.logical_and(x>xi, np.logical_and(x>xi*3.0, x<xi*(1+zem))))[0]
            tau_ls_dla[s_4] = tau_ls_dla[s_4] + a_dla_1 * (x[s_4]/xi)**2.
            tau_ls_dla[s_5] = tau_ls_dla[s_5] + a_dla_2 * (x[s_5]/xi)**3.

        # Eq. 25
        if zem < 1.2:
            s_6 = np.where(np.logical_and(x>x_ll, x<x_ll*(1+zem)))[0]
            tau_lc_laf[s_6] = 0.325*((x[s_6]/x_ll)**1.2 - (1+zem)**-0.9 * (x[s_6]/x_ll)**2.1)

        # Eq. 26
        elif zem < 4.7:
            s_7 = np.where(np.logical_and(x>x_ll, np.logical_and(x<2.2*x_ll, x<x_ll*(1+zem))))[0]
            s_8 = np.where(np.logical_and(x>x_ll, np.logical_and(x>2.2*x_ll, x<x_ll*(1+zem))))[0]
            tau_lc_laf[s_7] = 2.55e-2*(1+zem)**1.6 * (x[s_7]/x_ll)**2.1 + 0.325*(x[s_7]/x_ll)**1.2 - 0.250*(x[s_7]/x_ll)**2.1
            tau_lc_laf[s_8] = 2.55e-2*((1+zem)**1.6 * (x[s_8]/x_ll)**2.1 - (x[s_8]/x_ll)**3.7)

        # Eq. 27
        else:
            s_9 = np.where(np.logical_and(x>x_ll, np.logical_and(x<2.2*x_ll, x<x_ll*(1+zem))))[0]
            s_a = np.where(np.logical_and(x>x_ll, np.logical_and(x>2.2*x_ll, np.logical_and(x<5.7*x_ll, x<x_ll*(1+zem)))))[0]
            s_b = np.where(np.logical_and(x>x_ll, np.logical_and(x>2.2*x_ll, np.logical_and(x>5.7*x_ll, x<x_ll*(1+zem)))))[0]
            tau_lc_laf[s_9] = 5.22e-4*(1+zem)**3.4 * (x[s_9]/x_ll)**2.1 + 0.325*(x[s_9]/x_ll)**1.2 - 3.14e-2*(x[s_9]/x_ll)**2.1
            tau_lc_laf[s_a] = 5.22e-4*(1+zem)**3.4 * (x[s_a]/x_ll)**2.1 + 0.218*(x[s_a]/x_ll)**2.1 - 2.55e-2*(x[s_a]/x_ll)**3.7
            tau_lc_laf[s_b] = 5.22e-2*((1+zem)**3.4 * (x[s_b]/x_ll)**2.1 - (x[s_b]/x_ll)**5.5)

        # Eq. 28
        if zem < 2.0:
            s_c = np.where(np.logical_and(x>x_ll, x<x_ll*(1+zem)))[0]
            tau_lc_dla[s_c] = 0.211*(1+zem)**2. - 7.66e-2*(1+zem)**2.3 * (x[s_c]/x_ll)**(-0.3) - 0.135*(x[s_c]/x_ll)**2.

        # Eq. 29
        else:
            s_d = np.where(np.logical_and(x>x_ll, np.logical_and(x<3.0**x_ll, x<x_ll*(1+zem))))[0]
            s_e = np.where(np.logical_and(x>x_ll, np.logical_and(x>3.0*x_ll, x<x_ll*(1+zem))))[0]
            tau_lc_dla[s_d] = 0.634 + 4.7e-2*(1+zem)**3. - 1.78e-2*(1+zem)**3.3 * (x[s_d]/x_ll)**(-0.3) \
                              - 0.135*(x[s_d]/x_ll)**2. - 0.291*(x[s_d]/x_ll)**(-0.3)

            tau_lc_dla[s_e] = 4.7e-2*(1+zem)**3. - 1.78e-2*(1+zem)**3.3 * (x[s_e]/x_ll)**(-0.3) - 2.92e-2*(x[s_e]/x_ll)**3.

        tau = tau_ls_laf + tau_ls_dla + tau_lc_laf + tau_lc_dla

        corr = np.ones(len(x))
        z = (x/xem_d['Ly_a']).value - 1
        corr[z<zem] = np.exp(tau[z<zem])

        if apply:
            taucorr = dc(self._t[input_col])
            taucorr *= corr
            self._t[input_col+'_taucorr'] = taucorr

        return corr

    def _region_extract(self, xmin, xmax, verbose=True):
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

        if len(reg.t) == 0 and verbose:
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

    def _x_convert(self, zem=0, xunit=au.km/au.s, _update_zem=True):
        if _update_zem:
            print('upd')
            self._zem = zem
        self._xunit = xunit
        self._xunit_old = self._t['x'].unit
        self.x = x_convert(self.x, zem, xunit)
        self.xmin = x_convert(self.xmin, zem, xunit)
        self.xmax = x_convert(self.xmax, zem, xunit)
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
        self.y = self.y * fact * self.y.unit
        self.dy = self.dy * fact * self.dy.unit
        return 0
