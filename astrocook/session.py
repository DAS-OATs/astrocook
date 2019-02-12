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
        self.seq = ['spec', 'lines']

    def convert_x(self, zem=0, xunit=au.km/au.s):
        """ @brief Convert the x axis to wavelength or velocity units.
        @param zem Emission redshift, to use as a 0-point for velocities
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        try:
            zem = float(zem)
        except:
            print(prefix, msg_param_fail)
        xunit = au.Unit(xunit)

        for s in self.seq:
            try:
                getattr(self, s)._convert_x(zem, xunit)
            except:
                pass
        return 0

    def convert_y(self, yunit=au.electron/au.nm):
        """ @brief Convert the x axis to wavelength or velocity units.
        @param zem Emission redshift, to use as a 0-point for velocities
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        yunit = au.Unit(yunit)

        for s in self.seq:
            try:
                getattr(self, s)._convert_y(yunit=yunit)
            except:
                pass
        return 0

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
        kwargs = {'path': path, 'name': name}
        for s in self.seq:
            try:
                kwargs[s] = getattr(self, s)._extract_region(xmin, xmax)
            except:
                kwargs[s] = None
        if kwargs['spec'] != None:
            new = Session(**kwargs)
        else:
            new = None
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

        # DRS spectrum
        if instr == 'ESPRESSO' and catg[0:3] == 'S1D':
            self.spec = format.espresso_drs_spectrum(hdul)

        # DAS spectrum
        if instr == 'ESPRESSO' and catg[1:5] == 'SPEC':
            self.spec = format.espresso_das_spectrum(hdul)

        # ESO-MIDAS spectrum
        if orig == 'ESO-MIDAS':
            self.spec = format.eso_midas(hdul)
