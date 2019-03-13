from .format import Format
from .message import *
from .model import Model
from .line_list import LineList
from .spectrum import Spectrum
from .syst_list import SystList
from .vars import *
from astropy import units as au
from astropy.io import ascii, fits
from copy import deepcopy as dc
import numpy as np

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

    def _append(self, frame, append=True):
        if append and hasattr(self, frame.__name__):
            getattr(self, frame.__name__)._append(frame)
        else:
            setattr(self, frame.__name__, frame)


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

    def convolve_gauss(self, std=20, input_col='y', output_col='conv'):
        """@brief Convolve a spectrum column with a profile using FFT transform.
        @param std Standard deviation of the gaussian (km/s)
        @param input_col Input column
        @param output_col Output column
        @return 0
        """

        try:
            std = float(std) * au.km/au.s
        except:
            print(prefix, msg_param_fail)

        self.spec._convolve_gauss(std, input_col, output_col)
        return 0

    def extract_nodes(self, delta_x=1500, xunit=au.km/au.s):
        """ @brief Extract nodes from a spectrum. Nodes are averages of x and y
        in slices, computed after masking lines.
        @param delta_x Size of slices
        @param xunit Unit of wavelength or velocity
        @return 0
        """
        try:
            xunit = au.Unit(xunit)
            delta_x = float(delta_x)*xunit
        except:
            print(prefix, msg_param_fail)

        x, xmin, xmax, y, dy = self.spec._extract_nodes(delta_x, xunit)
        self.nodes = Spectrum(x, xmin, xmax, y, dy, self.spec._xunit,
                              self.spec._yunit)

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

        kwargs = {'path': self.path, 'name': self.name}
        for s in self.seq:
            try:
                #kwargs[s] = dc(getattr(self, s)._extract_region(xmin, xmax))
                kwargs[s] = getattr(self, s)._extract_region(xmin, xmax)
            except:
                try: # For Model
                    #kwargs[s] = dc(getattr(self, s)) # For Model
                    kwargs[s] = getattr(self, s) # For Model
                except:
                    kwargs[s] = None
        if kwargs['spec'] != None:
            new = Session(**kwargs)
        else:
            new = None

        return new

    def interp_nodes(self, smooth=0):
        """ @brief Interpolate nodes with a univariate spline to estimate the
        emission level.
        @param smooth Smoothing of the spline
        @return 0
        """
        if not hasattr(self, 'nodes'):
            print(prefix, "I need nodes to interpolate. Please try Spectrum > "
                  "Extract nodes first.")
            return None

        smooth = float(smooth)

        self.spec._interp_nodes(self.lines, self.nodes)
        return 0


    def find_peaks(self, col='conv', kind='min', kappa=3.0, append=True):
        """ @brief Find the peaks in a spectrum column. Peaks are the extrema
        (minima or maxima) that are more prominent than a given number of
        standard deviations. They are saved as a list of lines.
        @param col Column where to look for peaks
        @param kind Kind of extrema ('min' or 'max')
        @param kappa Number of standard deviations
        @param append Append peaks to existing line list
        @return 0
        """

        spec = self.spec
        if col not in spec.t.colnames:
            print(prefix, "The spectrum has not a column named '%s'. Please "\
                  "pick another one." % col)
            return None
        kappa = float(kappa)

        peaks = spec._find_peaks(col, kind, kappa)

        lines = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.y, peaks.dy,
                         spec._xunit, spec._yunit, spec._meta)

        if append and self.lines != None:
            self.lines._append(lines)
        else:
            self.lines = lines

        self.lines_kind = 'peaks'
        spec._mask_lines(self.lines)
        print(prefix, "I'm using peaks as lines.")

        return 0

    def add_fit(self, series='CIV', z=1.6971, logN=13, b=10, resol=45000,
                append=True):
        """ @brief Add and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param append Append system to existing system list
        @return 0
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        systs = SystList(sess=self)
        if append and self.systs != None:
            self.systs._append(systs)
        else:
            self.systs = systs

        self.systs._add_fit(series, z, logN, b, resol)

        return 0

    def add_fit_from_lines(self, series='CIV', z_start=1.71, z_end=1.18,
                           dz=1e-4, logN=14, b=10, append=True):
        """ @brief Add and fit Voigt models to a line list, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param N Guess column density
        @param b Guess doppler broadening
        @param append Append systems to existing system list
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        dz = float(dz)
        logN = float(logN)
        b = float(b)

        z_range = self.lines._syst_cand(series, z_start, z_end, dz)

        systs = SystList(sess=self)
        systs._add_fit(series, z_range, logN, b, 45000, verb=True)
        if append and self.systs != None:
            #print(self.systs._t)
            #print(systs)
            self.systs._append(systs)
            #print(self.systs._t)
        else:
            self.systs = systs

        return 0

    def test_fit_slide(self, series='CIV', z_start=1.13, z_end=1.71, z_step=5e-4,
                       logN=14, b=10.0, append=True):
        """ @brief Slide a set of Voigt models across a spectrum and fit them
        where they suit the spectrum.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift step
        @param logN Column density (logarithmic)
        @param b Doppler parameter
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        logN = float(logN)
        b = float(b)

        z_range = np.arange(z_start, z_end, z_step)

        systs = SystList(sess=self)
        systs._test_fit(self.spec, series, z_range, logN, b, 45000, verb=True)
        """
        if append and self.systs != None:
            #print(self.systs._t)
            #print(systs)
            self.systs._append(systs)
            #print(self.systs._t)
        else:
            self.systs = systs
        """
        return 0

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
