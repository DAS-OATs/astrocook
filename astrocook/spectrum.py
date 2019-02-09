from .frame import Frame
from .message import *
from astropy import units as au
from astropy import constants as aconst
from astropy import table as at
from copy import deepcopy as dc
import numpy as np
#from scipy.signal import fftconvolve

prefix = "Spectrum:"

class Spectrum(Frame):
    """Class for spectra

    A Spectrum is a Frame with methods for handling spectral operations."""

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
        super(Spectrum, self).__init__(x, xmin, xmax, y, dy, xunit, yunit, meta,
                                       dtype)

    def convolve_gauss(self, std=20):
        """@brief Convolve a spectrum with a profile using FFT transform and
        produces a new spectrum

        @param std Standard deviation of the gaussian
        @return Convolved spectrum.
        """

        # Create profile
        self.convert_x()
        x = self.x
        mean = np.median(x)*x.unit
        prof = np.exp(-((x - mean) / std).value**2)
        if (len(prof) % 2 == 0):
            prof = prof[:-1]
        prof = prof / np.sum(prof)

        # Convolve
        yconv = dc(self.y)
        #yconv = self._safe(yconv)
        self._t['conv'] = yconv# fftconvolve(yconv, prof, mode=mode)*yunit
        print(self._t)


    def convert_x(self, zem=0, xunit=au.km/au.s):
        """@brief Convert wavelengths into velocities and vice versa
        @param zem Emission redshift
        @param xunit Unit of velocity or wavelength
        """

        try:
            zem = float(zem)
        except:
            print(prefix, msg_param_fail)
        xunit = au.Unit(xunit)

        xem = (1+zem) * 121.567*au.nm
        equiv = [(au.nm, au.km / au.s,
                  lambda x: np.log(x/xem.value)*aconst.c.to(au.km/au.s),
                  lambda x: np.exp(x/aconst.c.to(au.km/au.s).value)*xem.value)]

        self._xunit = xunit
        self.x = self.x.to(xunit, equivalencies=equiv)
        self.xmin = self.xmin.to(xunit, equivalencies=equiv)
        self.xmax = self.xmax.to(xunit, equivalencies=equiv)
        return 0

    def extract_region(self, xmin, xmax):
        """ @brief Extract a spectral region as a new spectrum
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return Spectral region
        """

        try:
            xmin = float(xmin) * au.nm
            xmax = float(xmax) * au.nm
        except:
            print(prefix, msg_param_fail)
            return None

        if xmin > xmax:
            temp = xmin
            xmin = xmax
            xmax = temp
            print(prefix, msg_param_swap)


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
