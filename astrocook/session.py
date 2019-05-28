from . import version
from .cookbook import Cookbook
from .format import Format
from .functions import detect_local_minima
from .line_list import LineList
from .message import *
from .model import Model
from .spectrum import Spectrum
from .syst_list import SystList
from .syst_model import SystModel
#from .model_list import ModelList
from .vars import *
from astropy import units as au
from astropy.io import ascii, fits
from copy import deepcopy as dc
from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.signal import argrelmin
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
        self.cb = Cookbook(self)

    def _append(self, frame, append=True):
        if append and hasattr(self, frame.__name__):
            getattr(self, frame.__name__)._append(frame)
        else:
            setattr(self, frame.__name__, frame)

    """
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
    """

    def add_syst(self, series='Ly_a', z=2.0, logN=14, b=10, resol=70000,
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

        self.cb._append_syst()
        self.cb._fit_syst(series, z, logN, b, resol, maxfev)
        self.systs._clean(chi2r_thres)
        self.cb._update_spec()

        return 0


    def add_syst_from_lines(self, series='Ly_a', z_start=0, z_end=6,
                            dz=1e-4, logN=13, b=10, resol=45000,
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

        self.cb._append_syst()

        for i, z in enumerate(z_range):
            self.cb._fit_syst(series, z, logN, b, resol, maxfev)
            print(prefix, "I've fitted a %s system at redshift %2.4f (%i/%i)..."\
                  % (series, z, i+1, len(z_range)), end='\r')
        print(prefix, "I've fitted %i %s systems between redshift %2.4f and "\
              "%2.4f." % (len(z_range), series, z_range[0], z_range[-1]))

        self.systs._clean(chi2r_thres)
        self.cb._update_spec()

        return 0


    def add_syst_from_resids(self, z_start=1.71, z_end=1.18, dz=5e-5,
                             resol=45000, logN=13, b=5, chi2r_thres=2.0,
                             maxfev=100):
        """ @brief Add and fit Voigt models from residuals of previously
        fitted models.
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param N Guess column density
        @param b Guess doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to find models to improve
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

        systs = self.systs

        old = systs._t[np.where(systs._t['chi2r'] > chi2r_thres)]

        for i, o in enumerate(old):
            o_z = o['z']
            o_series = o['series']
            o_id = o['id']
            zmin = o_z-1e-3
            zmax = o_z+1e-3
            xmin = (1.+zmin)*xem_d[series_d[o_series][0]]
            xmax = (1.+zmax)*xem_d[series_d[o_series][-1]]


            chi2r_old = np.inf
            count = 0
            while True:

                spec = dc(self.spec)
                spec._convolve_gauss(std=10, input_col='deabs', verb=False)
                reg = spec._extract_region(xmin, xmax)
                peaks = reg._find_peaks()
                resids = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.y,
                                  peaks.dy, reg._xunit, reg._yunit, reg._meta)
                z_sel = resids._syst_cand(o_series, z_start, z_end, dz)

                # If no residuals are found, add a system at the init. redshift
                if len(z_sel)==0:
                    z_cand = o_z
                else:
                    z_cand = z_sel[0]

                print(prefix, "I'm improving a %s system at redshift %2.4f "\
                      "(%i/%i): trying to add a component at redshift %2.4f..."\
                      % (o_series, o_z, i+1, len(old), z_cand), end='\r')

                t_old, mods_t_old = systs._freeze()

                self._mods_t_old = mods_t_old
                #systs._append(SystList(id_start=len(self.systs._t)),
                systs._append(SystList(id_start=np.max(self.systs._t['id'])+1),
                              unique=False)
                self.cb._fit_syst(o_series, z_cand, logN, b, resol, maxfev)

                if systs._t['chi2r'][systs._t['id']==o_id]>chi2r_old:
                    systs._unfreeze(t_old, mods_t_old)
                    break

                chi2r_old = systs._t['chi2r'][systs._t['id']==o_id]
                self.cb._update_spec()
                count += 1

                if systs._t['chi2r'][systs._t['id']==o_id]<chi2r_thres: break

            if count == 0:
                print(prefix, "I've not improved the %s system at redshift "\
                      "%2.4f (%i/%i): I was unable to add useful components."\
                      % (o_series, o_z, i+1, len(old)))
            else:
                print(prefix, "I've improved a %s system at redshift %2.4f "\
                      "(%i/%i) by adding %i components.                       "\
                      "  " % (o_series, o_z, i+1, len(old), count))


        return 0


    def add_syst_slide(self, series='CIV',
                       z_start=0, z_end=6, z_step=5e-4,
                       logN_start=12, logN_end=11, logN_step=-0.1,
                       b_start=2, b_end=22, b_step=2,
                       resol=45000, col='y', chi2r_thres=2, maxfev=100):
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

        # Create x-swapped spectrum to monitor correctness
        sess = dc(self)
        sess.spec._t['x'] = sess.spec._t['x'][::-1]

        # Previously fitted systems are left fixed...
        systs_old = dc(self.systs)

        if hasattr(systs_old, '_t'):
            #self.systs = SystList(id_start=len(systs_old._t))
            self.systs = SystList(id_start=np.max(systs_old._t['id'])+1)
        else:
            self.systs = SystList()
        chi2a = np.full((len(logN_range),len(b_range),len(z_range)), np.inf)
        self.corr = np.empty((len(logN_range),len(b_range)))
        self.corr_logN = logN_range
        self.corr_b = b_range
        self.corr_e = (b_range[0]-b_step*0.5, b_range[-1]+b_step*0.5,
                       logN_range[0]-logN_step*0.5,logN_range[-1]+logN_step*0.5)
        for ilogN, logN in enumerate(logN_range):
            for ib, b in enumerate(b_range):
                xm, ym, ym_0, ym_1, ym_2 = self.cb._create_doubl(series, z_mean,
                                                                 logN, b, resol)
                cond_c = 0
                cond_swap_c = 0
                for iz, z in enumerate(z_range):
                    print(prefix, "I'm testing a %s system (logN=%2.2f, "
                          "b=%2.2f) at redshift %2.4f (%i/%i)..." \
                          % (series, logN, b, z, iz+1, len(z_range)), end='\r')
                    self.spec._shift_rf(z)
                    sess.spec._shift_rf(z)
                    #systs = SystList()
                    cond, chi2, chi2_0 = \
                        self.cb._test_doubl(xm, ym, ym_0, ym_1, ym_2, col)
                    cond_swap, _, _ = \
                        sess.cb._test_doubl(xm, ym, ym_0, ym_1, ym_2, col)
                    if cond:
                        chi2a[ilogN, ib, iz] = chi2
                        cond_c += 1
                    if cond_swap:
                        cond_swap_c += 1
                self.corr[ilogN, ib] = \
                    max(1-np.array(cond_swap_c)/np.array(cond_c), 0)
                print(prefix, "I've tested a %s system (logN=%2.2f, "\
                      "b=%2.2f) between redshift %2.4f and %2.4f and found %i "\
                      "coincidences"
                      % (series, logN, b, z_range[0], z_range[-1], cond_c),
                      end='')
                if cond_c != 0:
                    #print(cond_c, cond_swap_c)
                    print(" (estimated correctness=%2.0f%%)."
                          % (100*self.corr[ilogN, ib]))
                else:
                    print(".")
        self.spec._shift_rf(0)

        # Find candidates choosing the local minima of chi2r at coincidences
        lm = np.logical_and(chi2a < np.inf, detect_local_minima(chi2a))
        chi2m = np.where(lm)

        print(prefix, "I've selected %i candidates among the coincidences."\
              % len(chi2m[0]))
        if maxfev > 0:
            for i in range(len(chi2m[0])):
                z = z_range[chi2m[2][i]]
                logN = logN_range[chi2m[0][i]]
                b = b_range[chi2m[1][i]]
                self.cb._fit_syst(series, z, logN, b, resol, maxfev)
                print(prefix, "I've fitted a %s system at redshift %2.4f "
                      "(%i/%i)..." % (series, z, i+1, len(chi2m[0])), end='\r')
            if len(chi2m[0]) > 0:
                print(prefix, "I've fitted %i %s systems between redshift "
                      "%2.4f and %2.4f."
                      % (len(chi2m[0]), series, z_range[chi2m[2][0]],
                         z_range[chi2m[2][-1]]))
                self.systs._clean(chi2r_thres)

        # ...and then appended
        if hasattr(systs_old, '_t'):
            if len(self.systs._t) > 0:
                self.systs._append(systs_old)
            else:
                self.systs = systs_old

        self.cb._update_spec()

        # Save the correctness as a two-entry table - to be modified
        self.corr_save = np.concatenate((np.array([logN_range]).T, self.corr),
                                         axis=-1)
        self.corr_save = np.concatenate((np.array([np.append(['nan'], b_range)]),
                                         self.corr_save), axis=0)

        return 0


    def compl_syst(self, series='CIV', n=100,
                   z_start=0, z_end=6, z_step=5e-4,
                   logN_start=12, logN_end=11, logN_step=-0.1,
                   b_start=2, b_end=22, b_step=2,
                   resol=45000, col='y', chi2r_thres=2, maxfev=100):
        """ @brief Estimate the completeness of system detection by simulating
        systems at random redshifts and sliding Voigt models to fit them
        @param series Series of transitions
        @param n Number of simulated realizations
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
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param maxfev Maximum number of function evaluation
        @return 0
        """

        n = int(n)
        z_start = float(z_start)
        z_end = float(z_end)
        logN_start = float(logN_start)
        logN_end = float(logN_end)
        logN_step = float(logN_step)
        b_start = float(b_start)
        b_end = float(b_end)
        b_step = float(b_step)
        resol = float(resol)
        chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        logN_range = np.arange(logN_start, logN_end, logN_step)
        b_range = np.arange(b_start, b_end, b_step)

        # Previously fitted systems are left fixed...
        sess = dc(self)

        z_min = np.min([(np.min(self.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_max = np.max([(np.max(self.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_mean = 0.5*(z_min+z_max)
        self.compl = np.empty((len(logN_range),len(b_range)))
        self.compl_e = (b_range[0]-b_step*0.5, b_range[-1]+b_step*0.5,
                        logN_range[0]-logN_step*0.5,
                        logN_range[-1]+logN_step*0.5)
        import cProfile
        import pstats
        pr = cProfile.Profile()
        pr.enable()
        for ilogN, logN in enumerate(logN_range):
            for ib, b in enumerate(b_range):

                cond_c = 0
                n_ok = 0
                #for r in range(n):
                sess.systs = dc(self.systs)
                sess.cb._append_syst()
                sess.systs._add(series, z_mean, logN, b, resol)
                xm, ym, _, _, _ = sess.cb._create_doubl(series, z_mean, logN, b,
                                                        resol)
                while n_ok < n:
                    print(prefix, "I'm estimating completeness of %s system "
                          "(logN=%2.2f, b=%2.2f, realization %i/%i)..."
                          #% (series, logN, b, r+1, n), end='\r')
                          % (series, logN, b, n_ok+1, n), end='\r')
                    z_rand = np.random.rand()*(z_end-z_start)+z_start
                    sess.spec = dc(self.spec)
                    """
                    sess.systs = dc(self.systs)
                    sess.cb._append_syst()
                    fail = sess.cb._simul_syst(series, z_rand, logN, b, resol,
                                               col)
                    """
                    sess.spec._shift_rf(z_rand)
                    fail = sess.cb._apply_doubl(xm, ym)
                    if not fail:
                        n_ok += 1
                        z_round = round(z_rand, 4)
                        xm, ym, ym_0, ym_1, ym_2 = sess.cb._create_doubl(
                            series, z_mean, logN, b, resol)
                        for iz, z in enumerate(np.arange(
                            z_round-1*z_step, z_round+1*z_step, z_step)):
                            sess.spec._shift_rf(z)
                            systs = SystList()
                            cond, chi2, chi2_0 = sess.cb._test_doubl(
                                xm, ym, ym_0, ym_1, ym_2, col)
                            if cond and np.abs(z-z_rand)<z_step:
                                cond_c += 1
                        sess.spec._shift_rf(0)

                self.compl[ilogN, ib] = cond_c/n_ok
                print(prefix, "I've estimated completeness of %s system "
                      "(logN=%2.2f, b=%2.2f) as %2.0f%%.                       "
                      % (series, logN, b, 100*self.compl[ilogN, ib]))

        pr.disable()
        ps = pstats.Stats(pr)
        ps.sort_stats('cumulative').print_stats(10)
        #self = dc(sess_old)
        # Save the completeness as a two-entry table - to be modified
        self.compl_save = np.concatenate((np.array([logN_range]).T, self.compl),
                                         axis=-1)
        self.compl_save = np.concatenate((np.array([np.append(['nan'], b_range)]),
                                          self.compl_save), axis=0)
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


    def open(self):

        format = Format()
        if self.path[-3:] == 'acs':
            root = '/'.join(self.path.split('/')[:-1])
            with tarfile.open(self.path) as arch:
                arch.extractall(path=root)
                hdul = fits.open(self.path[:-4]+'_spec.fits')
                hdr = hdul[1].header
        else:
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
            if 'XSHOOTER_REDUCE' in hist:
                instr = 'XSHOOTER'
                orig = 'REDUCE'
        except:
            pass

        if instr == None:
            print(prefix, "INSTRUME not defined.")
        if catg == None:
            print(prefix, "HIERARCH ESO PRO CATG not defined.")
        if orig == None:
            print(prefix, "ORIGIN not defined.")

        # Astrocook structures
        if orig == 'Astrocook':
            for s in self.seq:
                try:
                    hdul = fits.open(self.path[:-4]+'_'+s+'.fits')
                    setattr(self, s, format.astrocook(hdul, s))
                except:
                    pass

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

        # XSHOOTER_REDUCE spectrum
        if instr == 'XSHOOTER' and orig == 'REDUCE':
            hdul_e = fits.open(self.path[:-5]+'e.fits')
            self.spec = format.xshooter_reduce_spectrum(hdul, hdul_e)


    def save(self, path):

        root = path[:-4]
        stem = root.split('/')[-1]
        with tarfile.open(root+'.acs', 'w:gz') as arch:
            for s in self.seq:
                try:
                    if s=='systs':
                        np.savetxt(root+'_compl.dat', self.compl_save, fmt='%s')
                        np.savetxt(root+'_corr.dat', self.corr_save, fmt='%s')
                    name = root+'_'+s+'.fits'
                    obj = dc(getattr(self, s))
                    t = obj._t
                    t.meta = obj._meta
                    t.meta['ORIGIN'] = 'Astrocook'
                    t.meta['HIERARCH ASTROCOOK VERSION'] = version
                    t.meta['HIERARCH ASTROCOOK STRUCT'] = s
                    for c in t.colnames:
                        t[c].unit = au.dimensionless_unscaled
                    t.write(name, format='fits', overwrite=True)
                    arch.add(name, arcname=stem+'_'+s+'.fits')
                    os.remove(name)

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

    def simul_syst(self, series='Ly_a', z=2.0, logN=14, b=10, resol=70000,
                   col='y'):
        """ @brief Simulate a system by adding Voigt model onto a spectrum.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @return 0
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        self.cb._append_syst()
        self.cb._simul_syst(series, z, logN, b, resol, col)

        return 0
