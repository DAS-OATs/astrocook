from .format import Format
from .message import *
from .model import Model
from .line_list import LineList
from .spectrum import Spectrum
from .syst_list import SystList, SystList2
from .syst_model import SystModel, SystModel2
#from .model_list import ModelList
from .vars import *
from astropy import units as au
from astropy.io import ascii, fits
from copy import deepcopy as dc
import numpy as np
import tarfile

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


    def _fit_syst(self, series='CIV', z=2, logN=13, b=10, resol=70000,
                  maxfev=100, verb=True):

        self.systs._add2(series, z, logN, b, resol)
        mod = SystModel2(self.spec, self.systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)
        mod._fit(fit_kws={'maxfev': maxfev}, update=False)
        self.systs._update(mod)

        return 0


    def _test_doubl(self, series='Ly_a', z_mean=2.0, logN=14, b=10,
                    resol=70000):

        spec = self.spec
        spec._shift_rf(z_mean)
        mod = SystModel2(self.spec, self.systs, z0=0)
        mod._new_voigt(series, 0, logN, b, resol)
        spec._shift_rf(0.0)
        xm = mod._xf
        hlenm = len(xm)//2
        ym = mod.eval(x=xm, params=mod._pars)
        ym_0 = np.ones(len(xm))
        ym_1 = np.concatenate([ym[:-hlenm], np.ones(hlenm)])
        ym_2 = np.concatenate([np.ones(hlenm), ym[hlenm:]])

        return xm, ym, ym_0, ym_1, ym_2

    def _update_spec(self):

        spec = self.spec

        self.systs._xs = np.array(spec._safe(spec.x).to(au.nm))
        s = spec._where_safe

        y = spec.y
        if 'model' not in spec._t.colnames:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*y.unit
        if 'deabs' not in spec._t.colnames:
            spec._t['deabs'] = y

        #s = self.systs._s
        cont = spec._t['cont']
        model = spec._t['model']
        deabs = spec._t['deabs']

        model[s] = cont[s]
        for i, r in enumerate(self.systs._mods_t):
            mod = r['mod']
            model[s] = mod.eval(x=self.systs._xs, params=mod._pars) * model[s]
        deabs[s] = cont[s] + y[s] - model[s]


    def add_syst(self, series='CIV', z=1.6971, logN=13, b=10, resol=70000,
                 chi2r_thres=np.inf, maxfev=100):
        """ @brief Add and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param maxfev Maximum number of function evaluation
        @return 0
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)
        chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        systs = SystList2()
        if self.systs != None:
            self.systs._append(systs)
        else:
            self.systs = systs
        self._fit_syst(series, z, logN, b, resol, maxfev)
        self.systs._clean(chi2r_thres)
        print(prefix, "I've fitted a %s systems between at redshift %2.4f."\
              % (series, z))

        self._update_spec()

        return 0


    def add_syst_from_lines(self, series='CIV', z_start=1.71, z_end=1.18,
                            dz=1e-4, logN=14, b=10, resol=70000,
                            chi2r_thres=np.inf, maxfev=100):
        """ @brief Add and fit Voigt models to a line list, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param N Guess column density
        @param b Guess doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param maxfev Maximum number of function evaluation
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        dz = float(dz)
        logN = float(logN)
        b = float(b)
        chi2r_thres = float(chi2r_thres)
        resol = float(resol)
        maxfev = int(maxfev)

        z_range = self.lines._syst_cand(series, z_start, z_end, dz)

        systs = SystList2()
        if self.systs != None:
            self.systs._append(systs)
        else:
            self.systs = systs

        for i, z in enumerate(z_range):
            self._fit_syst(series, z, logN, b, resol, maxfev)
            print(prefix, "I've fitted a %s system at redshift %2.4f (%i/%i)…"\
                  % (series, z, i+1, len(z_range)), end='\r')
        print(prefix, "I've fitted %i %s systems between redshift %2.4f and "\
              "%2.4f." % (len(z_range), series, z_range[0], z_range[-1]))

        self.systs._clean(chi2r_thres)
        self._update_spec()

        return 0


    def add_syst_slide(self, series='CIV',
                       z_start=1.13, z_end=1.71, z_step=5e-4,
                       logN_start=14, logN_end=14.1, logN_step=0.1,
                       b_start=10, b_end=15, b_step=5,
                       resol=70000, col='deabs',
                       chi2_fact=1.0, chi2r_thres=2.0, maxfev=100):
        """ @brief Slide a set of Voigt models across a spectrum and fit them
        where they suit the spectrum.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift step
        @param logN_start Start column density (logarithmic)
        @param logN_end End column density (logarithmic)
        @param logN_step Column density step (logarithmic)
        @param b_start Start Doppler parameter
        @param b_end End Doppler parameter
        @param b_step Doppler parameter step
        @param resol Resolution
        @param col Column where to test the models
        @param chi2_fact Maximum ratio between the chi2 of model and null case
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param maxfev Maximum number of function evaluation
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        logN_start = float(logN_start)
        logN_end = float(logN_end)
        logN_step = float(logN_step)
        b_start = float(b_start)
        b_end = float(b_end)
        b_step = float(b_step)
        resol = float(resol)
        chi2_fact = float(chi2_fact)
        chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        z_range = np.arange(z_start, z_end, z_step)
        z_min = np.min([(np.min(self.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_max = np.max([(np.max(self.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_range = z_range[np.where(np.logical_and(z_range > z_min,
                                                  z_range < z_max))]
        z_mean = 0.5*(z_min+z_max)
        logN_range = np.arange(logN_start, logN_end, logN_step)
        b_range = np.arange(b_start, b_end, b_step)

        # Previously fitted systems are left fixed...
        systs_old = dc(self.systs)

        self.systs = SystList2()
        for logN in logN_range:
            for b in b_range:
                xm, ym, ym_0, ym_1, ym_2 = self._test_doubl(series, z_mean, logN, b, resol)
                z_true = []
                for z in z_range:
                    self.spec._shift_rf(z)
                    systs = SystList2()
                    cond, chi2, chi2_0 = \
                        systs._test(self.spec, xm, ym, ym_0, ym_1, ym_2,
                                    col, chi2_fact)
                    if cond:
                        z_true.append(z)

        self.spec._shift_rf(0)

        #self.systs = SystList2()
        for z in z_true:
            self._fit_syst(series, z, logN_start, b_start, resol, maxfev)
            #print(self.systs._mods_t)
        self.systs._clean(chi2r_thres)

        # ...and then appended
        self.systs._append(systs_old)

        self._update_spec()


        return 0


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

    """
    def test_fit_slide(self, series='CIV', z_start=1.13, z_end=1.71,
                       z_step=5e-4, logN=14, b=10.0, col='deabs', resol=70000,
                       chi2_fact=1.0, chi2r_thres=2.0, maxfev=100, append=True):

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        logN = float(logN)
        b = float(b)
        resol = float(resol)
        chi2_fact = float(chi2_fact)
        if chi2r_thres == None or chi2r_thres == 'None':
            chi2r_thres = np.inf
        else:
            chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        z_range = np.arange(z_start, z_end, z_step)

        systs = SystList(sess=self)
        #mods = ModelList()
        systs._test_fit(self.spec, series, z_range, logN, b, resol, col,
                        chi2_fact, chi2r_thres, fit_kws={'maxfev': maxfev},
                        verb=True)
        if append and self.systs != None:
            self.systs._append(systs)
            #self.systs._mods._append(systs._mods)
        else:
            self.systs = systs
        self._update_spec()

        return 0
    """

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

    def save(self, path):
        for s in self.seq:
            try:
                name = path[:-4]+'_'+s+'.fits'
                t = dc(getattr(self, s).t)
                for c in t.colnames:
                    t[c].unit = au.dimensionless_unscaled
                t.write(name, format='fits', overwrite=True)
            except:
                pass


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
