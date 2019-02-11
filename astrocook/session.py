from .format import Format
from .message import *
from .spectrum import Spectrum
from astropy import units as au
from astropy.io import fits

prefix = "Session:"

class Session(object):
    """ Class for sessions.

    A Session is a self-sufficient set of analysis operations."""

    def __init__(self,
                 path=None,
                 name=None,
                 spec=None,
                 lines=None):
        self.path = path
        self.name = name
        self.spec = spec
        self.lines = lines

    def extract_region(self, xmin, xmax):
        """ @brief Extract a spectral region as a new frame.
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

        path = self.path
        name = self.name
        spec = self.spec._extract_region(xmin, xmax)
        lines = self.lines._extract_region(xmin, xmax)
        new = Session(path, name, spec, lines)
        return new

    def open(self):
        format = Format()
        hdul = fits.open(self.path)
        try:
            instr = hdul[0].header['INSTRUME']
        except:
            instr = None
            print(prefix, "INSTRUME not defined.")
        try:
            catg = hdul[0].header['HIERARCH ESO PRO CATG']
        except:
            catg = None
            print(prefix, "HIERARCH ESO PRO CATG not defined.")
        try:
            orig = hdul[0].header['ORIGIN']
        except:
            orig = None
            print(prefix, "ORIGIN not defined.")

        # DAS spectrum
        if instr == 'ESPRESSO' and catg == 'RSPEC':
            self.spec = format.das_spectrum(hdul)
        if orig == 'ESO-MIDAS':
            self.spec = format.eso_midas(hdul)
