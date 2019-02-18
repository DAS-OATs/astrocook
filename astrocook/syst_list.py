from .model import Model
from .vars import *
from astropy import table as at
from astropy import units as au
import numpy as np

prefix = "System list:"

class SystList(object):
    """ Class for system lists

    A SystList is a list of absorption systems with methods for handling
    spectral lines. """

    def __init__(self,
                 sess,
                 series=[],
                 func=[],
                 z=[],
                 zmin=[],
                 zmax=[],
                 dtype=float):
        self._sess = sess
        t = at.Table()
        zunit = au.dimensionless_unscaled
        t['series'] = at.Column(np.array(func, ndmin=1), dtype='S100')
        t['func'] = at.Column(np.array(func, ndmin=1), dtype='S5')
        t['z'] = at.Column(np.array(z, ndmin=1), dtype=dtype, unit=zunit)
        t['zmin'] = at.Column(np.array(zmin, ndmin=1), dtype=dtype, unit=zunit)
        t['zmax'] = at.Column(np.array(zmax, ndmin=1), dtype=dtype, unit=zunit)
        self._t = t
        self._dtype = dtype

    @property
    def t(self):
        return self._t

    @property
    def series(self):
        return self._t['series']

    @property
    def func(self):
        return self._t['func']

    @property
    def z(self):
        return au.Quantity(self._t['z'])

    @property
    def zmin(self):
        return au.Quantity(self._t['zmin'])

    @property
    def zmax(self):
        return au.Quantity(self._t['zmax'])

    @series.setter
    def func(self, val):
        self._t['series'] = np.array(val, dtype='S100')

    @func.setter
    def func(self, val):
        self._t['func'] = np.array(val, dtype='S5')

    @z.setter
    def z(self, val, dtype=float):
        self._t['z'] = np.array(val, dtype=dtype)
        self._t['z'].unit = val.unit

    @zmin.setter
    def zmin(self, val, dtype=float):
        self._t['zmin'] = np.array(val, dtype=dtype)
        self._t['zmin'].unit = val.unit

    @zmax.setter
    def zmax(self, val, dtype=float):
        self._t['zmax'] = np.array(val, dtype=dtype)
        self._t['zmax'].unit = val.unit

    def fit_single(self, series='Ly_a', z=2.0, dz=0.005, N=1e14, b=20, btur=0,
                   resol=35000, cont_ampl=0.1):
        """ @brief Create a voigt model of a system
        @param series Series of transitions
        @param z Redshift
        @param dz Redshift tolerance
        @param N Column density
        @param b Doppler broadening
        @param btur Turbulence broadening
        @param resol Resolution
        @return 0
        """

        z = float(z)
        dz = float(dz)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)
        cont_ampl = float(cont_ampl)

        # Create model
        func = 'voigt'
        zmin = z-dz
        zmax = z+dz
        pars = {'N': N, 'b': b, 'b_tur': btur, 'resol': resol}
        self._t.add_row([series, func, z, zmin, zmax])#, params, vary, expr)
        model = Model(series=series, z=z, zmin=zmin, zmax=zmax, pars=pars)
        comp = model._comp

        # Fit model
        spec = self._sess.spec
        x = np.array(spec._safe(spec.x).to(au.nm))  # Full x array, without NaNs
        s = spec._where_safe
        c = np.array([], dtype=int)
        for t in series_d[series]:
            c = np.append(c, np.where(
                    np.logical_and(x > (1+zmin)*xem_d[t].value,
                                   x < (1+zmax)*xem_d[t].value))[0])
        xc = np.array(x[c])  # Intervals used for fitting
        if 'deabs' in spec._t.colnames:
            yc = np.array(spec._t['deabs'][c]/spec._t['cont'][c])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        wc = np.array(spec._t['cont'][c]/spec.dy[c])
        fit = comp.fit(yc, comp._pars, x=xc, weights=wc)

        # Update table
        if 'model' not in spec._t.colnames:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*spec.y.unit
            spec._t['model'][s] = spec._t['cont'][s]
        spec._t['model'][s] = fit.eval(x=x) * spec._t['model'][s]
        if 'deabs' not in spec._t.colnames:
            spec._t['deabs'] = spec.y
        spec._t['deabs'][s] = spec._t['cont'][s]+spec.y[s]-spec._t['model'][s]

        return 0
