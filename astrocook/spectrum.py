from .frame import Frame
from astropy import units as au
from copy import deepcopy as dc
import numpy as np

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

    def extract_region(self, xmin=0.0, xmax=0.0):
        """ @brief Extract a spectral region as a new spectrum
        @param xmin Minimum wavelength
        @param xmax Maximum wavelength
        """

        reg = dc(self)
        where = np.full(len(reg.x), True)
        s = np.where(np.logical_and(reg.x > xmin, reg.x < xmax))
        where[s] = False
        reg._t.remove_rows(where)
        return reg
