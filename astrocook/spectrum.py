from .frame import Frame
from .message import *
from astropy import units as au
from copy import deepcopy as dc
import numpy as np

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
