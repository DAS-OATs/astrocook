from .format import Format
from .message import *
from .model import Model
#from .spectrum import Spectrum
from astropy import units as au
from astropy.io import ascii, fits
from copy import deepcopy as dc

prefix = "Session:"

class Session(object):
    """ Class for sessions.

    A Session is a self-sufficient set of analysis operations."""

    def __init__(self,
                 path=None,
                 name=None,
                 spec=None,
                 spec_form=None,
                 nodes=None,
                 lines=None,
                 systs=None,
                 mods=None):
        self.path = path
        self.name = name
        self.spec = spec
        self.spec_form = spec_form
        self.nodes = nodes
        self.lines = lines
        self.systs = systs
        self.mods = mods
        self.seq = ['spec', 'nodes', 'lines', 'systs', 'mods']

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
                kwargs[s] = dc(getattr(self, s)._extract_region(xmin, xmax))
            except:
                try: # For Model
                    kwargs[s] = dc(getattr(self, s)) # For Model
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
        hdr = hdul[0].header
        try:
            instr = hdr['INSTRUME']
        except:
            instr = None
        try:
            orig = hdr['ORIGIN']
        except:
            orig = None
        try:
            catg = hdr['HIERARCH ESO PRO CATG']
        except:
            catg = None

        try:
            hist = [i.split(' ') for i in str(hdr['HISTORY']).split('\n')]
            hist = [i for j in hist for i in j]
            if 'UVES_popler:' in hist:
                instr = 'UVES'
                orig = 'POPLER'
        except:
            pass

        if instr == None:
            print(prefix, "INSTRUME not defined.")
        if catg == None:
            print(prefix, "HIERARCH ESO PRO CATG not defined.")
        if orig == None:
            print(prefix, "ORIGIN not defined.")

        # ESO-MIDAS spectrum
        if orig == 'ESO-MIDAS':
            self.spec = format.eso_midas(hdul)

        # ESPRESSO DRS spectrum
        if instr == 'ESPRESSO' and catg[0:3] == 'S1D':
            self.spec = format.espresso_drs_spectrum(hdul)
            self.spec_form = format.espresso_spectrum_format(
                ascii.read('espr_spec_form.dat'))

        # ESPRESSO DAS spectrum
        if instr == 'ESPRESSO' and catg[1:5] == 'SPEC':
            self.spec = format.espresso_das_spectrum(hdul)
            self.spec_form = format.espresso_spectrum_format(
                ascii.read('espr_spec_form.dat'))

        # UVES POPLER spectrum
        if instr == 'UVES' and orig == 'POPLER':
            self.spec = format.uves_popler_spectrum(hdul)

        #self.model = Model(self)

    def shift_from_rf(self, z=0):
        """ @brief Shift x axis to the original frame.
        @return 0
        """

        try:
            z = float(z)
        except:
            print(prefix, msg_param_fail)

        for s in self.seq:
            try:
                getattr(self, s)._shift_rf(z)
            except:
                pass
        return 0

    def shift_to_rf(self, z=0):
        """ @brief Shift x axis to the rest frame.
        @param z Redshift to use for shifting
        @return 0
        """

        try:
            z = float(z)
        except:
            print(prefix, msg_param_fail)

        for s in self.seq:
            try:
                getattr(self, s)._shift_rf(z)
            except:
                pass
        return 0
