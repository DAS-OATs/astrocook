import numpy as np
from astropy.io import fits as fits
from astropy import units as u
from specutils import extinction
import copy


class spec1d:
    """Class suitable to store a 1D spectrum"""

    def __init__(self, x, y, dy, 
                 dx=None, 
                 xUnit=u.dimensionless_unscaled, 
                 yUnit=u.dimensionless_unscaled, 
                 good=None,
                 exptime=float('nan'), order=-1, meta=None):
        ''' Constructor for the spec1d class. '''
        self._order   = int(order)
        self._exptime = float(exptime)
        if (meta is None):
            meta = {}
        self._meta    = copy.deepcopy(meta)
        self._x       = np.asarray(copy.deepcopy(x), dtype=float)
        if (dx is None):
            dx = [1]*len(x) #TODO: check this 
        self._dx      = np.asarray(copy.deepcopy(dx) , dtype=float)
        self._xUnit   = copy.deepcopy(xUnit)
        self._y       = np.asarray(copy.deepcopy(y) , dtype=float)
        self._dy      = np.asarray(copy.deepcopy(dy), dtype=float)
        self._yUnit   = copy.deepcopy(yUnit)
        if (good is None):
            good = [1]*len(x)
        self._good    = np.asarray(copy.deepcopy(good), dtype=int)

        self._useGood = False



    @property
    def useGood(self):
        """Tells whether x, y, etc.  getters return only data from channels flagged as good."""
        return self._useGood

    @useGood.setter
    def useGood(self, value):
        self._useGood = value
        if self._useGood:
            self._igood = np.argwhere(self._good > 0)

    @property
    def good(self):
        """Good flag for each spectrum channel."""
        return self._good

    def nchan(self):
        """Number of channels in the spectrum."""
        if self._useGood:
            return len(self._igood)
        else:
            return len(self._x)

    @property
    def order(self):
        """Spectrum order (-1 if not specified)."""
        return self._order

    @property
    def exptime(self):
        """Spectrum exposure time (in seconds)."""
        return self._exptime

    @property
    def meta(self):
        """Meta information for the spectrum."""
        return self._meta

    @property
    def x(self):
        """Quantities associated to spectrum channels (e.g. wavelength, frequency, energies, etc.)."""
        if self._useGood:
            return self._x[self._igood]
        else:
            return self._x

    @property
    def dx(self):
        """Widths of spectrum channels, in the same units as the spectrum channels."""
        if self._useGood:
            return self._dx[self._igood]
        else:
            return self._dx

    @property
    def xUnit(self):
        """Physical unit for the x property, to be expressed as an astropy unit."""
        return self._xUnit

    @property
    def y(self):
        """Quantities associated to spectrum intensities (e.g. flux density, luminosity density, nuF_nu, lambdaF_lambda, etc.)."""
        if self._useGood:
            return self._y[self._igood]
        else:
            return self._y

    @property
    def dy(self):
        """Uncertainties associated to y values."""
        if self._useGood:
            return self._dy[self._igood]
        else:
            return self._dy

    @property
    def yUnit(self):
        """Physical unit for the y property, to be expressed as an astropy unit."""
        return self._yUnit

    def convert(self, xUnit=None, yUnit=None):
        """Convert x and/or y values into equivalent quantities."""
        if not (xUnit is None):
            q = self._x * self._xUnit
            p = q.to(xUnit, equivalencies=u.spectral())
            self._x = p.value
            self._xUnit = xUnit
        
        if not (yUnit is None):
            for i in range(0, len(self._x)):
                q = self._y[i] * self._yUnit
                p = q.to(yUnit, equivalencies=u.spectral_density(self._x[i] * self._xUnit))
                self._y[i] = p.value

                q = self._dy[i] * self._yUnit
                p = q.to(yUnit, equivalencies=u.spectral_density(self._x[i] * self._xUnit))
                self._dy[i] = p.value
            self._yUnit = yUnit


    def deredden(self, A_v, model='od94'):
        extFactor = extinction.reddening(self._x * self._xUnit, A_v, model=model)
        self._y  *= extFactor
        self._dy *= extFactor
