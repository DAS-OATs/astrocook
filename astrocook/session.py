from . import version
from .cookbook import Cookbook
from .format import Format
from .functions import detect_local_minima
from .line_list import LineList
from .message import *
#from .model import Model
from .spectrum import Spectrum
from .syst_list import SystList
from .syst_model import SystModel
#from .model_list import ModelList
from .vars import *
#from astropy import constants as ac
from astropy import units as au
from astropy.io import ascii, fits
from copy import deepcopy as dc
from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.signal import argrelmin
import tarfile
import time

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
        """ @brief Fit a system
        @details Add and fit a Voigt model for a system.
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
        """ @brief Fit systems from line list
        @details Add and fit Voigt models to a line list, given a redshift
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
        if series == 'unknown':
            z_start = 0
            z_end = np.inf
        dz = float(dz)
        if logN is not None:
            logN = float(logN)
        b = float(b)
        chi2r_thres = float(chi2r_thres)
        resol = float(resol)
        maxfev = int(maxfev)

        if logN is None:
            z_range, logN_range = self.lines._syst_cand(series, z_start, z_end,
                                                        dz, logN=True)
            self.cb._append_syst()
            for i, (z, l) in enumerate(zip(z_range, logN_range)):
                ### Uncomment this
                self.cb._mod_syst(series, z, l, b, resol)

                ### Remove this
                #print(prefix, "I'm fitting a %s model at redshift %2.4f..."\
                #    % (series, z), end='\r')
                #self.cb._fit_syst(series, z, l, b, resol, maxfev)
                ###

            #print(z_range)
            #print(logN_range)
        else:
            z_range = self.lines._syst_cand(series, z_start, z_end, dz)
            self.cb._append_syst()
            #print(len(self.lines.t), len(z_range))

            #start = time.time()
            for i, z in enumerate(z_range):
                #print(i, z)
                ### Uncomment this
                self.cb._mod_syst(series, z, logN, b, resol)

                ### Remove this
                #print(prefix, "I'm fitting a %s model at redshift %2.4f..."\
                #    % (series, z), end='\r')
                #self.cb._fit_syst(series, z, logN, b, resol, maxfev)
                ###
            #end = time.time()

        ### Uncomment this
        #"""
        mods_t = self.systs._mods_t
        if len(z_range) > 0:
            print(prefix, "I've added %i %s system(s) in %i model(s) between "
                  "redshift %2.4f and %2.4f." % (len(z_range), series,
                  len(mods_t), z_range[0], z_range[-1]))

        if maxfev > 0:
            #print(mods_t['z0', 'id'])
            for i,m in enumerate(mods_t):
                print(prefix, "I'm fitting a %s model at redshift %2.4f "
                      "(%i/%i)..."\
                    % (series, m['z0'], i+1, len(mods_t)), end='\r')
                self.cb._fit_mod(m['mod'], maxfev)
            #print(self.systs._t)
            try:
                print(prefix, "I've fitted %i %s system(s) in %i model(s) "
                      "between redshift %2.4f and %2.4f." \
                      % (len(z_range), series, len(mods_t), z_range[0],
                         z_range[-1]))
            except:
                pass
        #"""
        if len(self.systs._t) > 0:
            self.systs._clean(chi2r_thres)
            self.cb._update_spec()
        else:
            self.systs = None

        return 0


    def add_syst_from_resids(self, z_start=0, z_end=6, dz=1e-4,
                             resol=45000, logN=11, b=5, chi2r_thres=1.0,
                             maxfev=100):
        """ @brief Fit systems from residuals
        @details Add and fit Voigt models from residuals of previously fitted
        models.
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param resol Resolution
        @param logN Guess column density
        @param b Guess doppler broadening
        @param chi2r_thres Reduced chi2 threshold to find models to improve
        @param maxfev Maximum number of function evaluation
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        dz = float(dz)
        resol = float(resol)
        logN = float(logN)
        b = float(b)
        chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        systs = self.systs

        #old = systs._t[np.where(systs._t['chi2r'] > chi2r_thres)]
        old = systs._mods_t[np.where(np.logical_or(
                systs._mods_t['chi2r'] > chi2r_thres,
                np.isnan(systs._mods_t['chi2r'])))]
        for i, o in enumerate(old):
            o_id = o['id'][0]
            o_series = systs._t[systs._t['id'] == o_id]['series'][0]
            o_z = np.array(systs._t['z'][systs._t['id']==o_id])[0]

            chi2r_old = np.inf
            #bic_old = np.inf
            count = 0
            count_good = 0

            while True:

                spec = dc(self.spec)
                spec._convolve_gauss(std=2, input_col='deabs', verb=False)
                try:
                #    ciao
                #except:
                #    print(systs._mods_t['mod'][i])
                    reg_x = systs._mods_t['mod'][i]._xm
                except:
                    break
                    #print("except")
                    #print(" ")
                    #reg_x = systs._mods_t['mod'][i-1]._xm
                    #reg_x = o['mod']._xm

                reg_xmin = np.interp(reg_x, spec.x.to(au.nm), spec.xmin.to(au.nm))
                reg_xmax = np.interp(reg_x, spec.x.to(au.nm), spec.xmax.to(au.nm))
                reg_y = np.interp(reg_x, spec.x.to(au.nm), spec.t['conv']-spec.t['cont'])
                #plt.scatter(reg_x, np.ones(len(reg_x)))
                reg_dy = np.interp(reg_x, spec.x.to(au.nm), spec.dy)
                reg = Spectrum(reg_x, reg_xmin, reg_xmax, reg_y, reg_dy)
                peaks = reg._find_peaks(col='y')#, mode='wrap')
                #print(peaks.t)

                resids = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.y,
                                  peaks.dy, reg._xunit, reg._yunit, reg._meta)
                #print(resids._t)
                #plt.scatter(peaks.x, peaks.y)
                #plt.show()
                z_cand = resids._syst_cand(o_series, z_start, z_end, dz,
                                           single=True)
                z_alt = resids._syst_cand('unknown', 0, np.inf, dz, single=True)
                #print(z_cand)
                #print(z_alt)


                # If no residuals are found, add a system at the init. redshift
                if z_cand == None:
                    #print("hey")
                    z_cand = o_z

                if z_alt == None:
                    #print("yeh")
                    z_alt = (1.+o_z)\
                            *xem_d[series_d[o_series][0]].to(au.nm).value

                # Randomize
                #z_cand = z_cand+np.random.normal(scale=0.0005)
                #z_alt = z_alt+np.random.normal(scale=0.0005)

                if count == 0:
                    t_old, mods_t_old = self.systs._freeze()
                systs._append(SystList(id_start=np.max(self.systs._t['id'])+1),
                              unique=False)
                cand = dc(self)
                alt = dc(self)

                cand_mod = cand.cb._fit_syst(o_series, z_cand, logN, b, resol, maxfev)
                alt_mod = alt.cb._fit_syst('unknown', z_alt, logN, b, resol, maxfev)
                """
                print(cand_mod._chi2r, alt_mod._chi2r,
                      cand_mod._aic, alt_mod._aic,
                      cand_mod._bic, alt_mod._bic)
                """
                self.systs._mods_t = dc(mods_t_old)
                #chi2r_cand = cand.systs._t['chi2r'][cand.systs._t['id']==o_id][0]
                #chi2r_alt = alt.systs._t['chi2r'][alt.systs._t['id']==o_id][0]
                chi2r_cand = cand_mod._chi2r
                chi2r_alt = alt_mod._chi2r
                """
                if chi2r_cand > chi2r_alt:#*2:#1.1:# and count > 3:
                    mod = self.cb._fit_syst('unknown', z_alt, logN, b, resol, maxfev)
                    #self.systs._add('unknown', z_alt, logN, b, resol)
                    #self.systs._update(alt_mod)
                    msg = "added an unknown component at wavelength %2.4f" \
                          % z_alt
                    chi2r = chi2r_alt
                else:
                    mod = self.cb._fit_syst(o_series, z_cand, logN, b, resol, maxfev)
                    #self.systs._add(o_series, z_cand, logN, b, resol)
                    #self.systs._update(cand_mod)
                    msg = "added a %s component at redshift %2.4f" \
                          % (o_series, z_cand)
                    chi2r = chi2r_cand
                """
                #print(chi2r_old, cand_mod._chi2r, alt_mod._chi2r)
                #if chi2r>=chi2r_old*100:# and count > 2:
                #if mod._chi2r>=chi2r_old:# and chi2r < 10:#50:
                #if ((cand_mod._chi2r>=chi2r_old and alt_mod._chi2r>= chi2r_old)
                #    and count_good > 5) \
                #    or (cand_mod._chi2r<10 or alt_mod._chi2r<10):
                if cand_mod._chi2r>=chi2r_old and alt_mod._chi2r>= chi2r_old:
                    self.systs._unfreeze(t_old, mods_t_old)
                    #mod._bic = bic_old
                    count += 1
                    #if chi2r < 10:
                    break
                else:
                    t_old, mods_t_old = self.systs._freeze()
                    if chi2r_cand > chi2r_alt:#*2:#1.1:# and count > 3:
                        mod = self.cb._fit_syst('unknown', z_alt, logN, b, resol, maxfev)
                        #self.systs._add('unknown', z_alt, logN, b, resol)
                        #self.systs._update(alt_mod)
                        msg = "added an unknown component at wavelength %2.4f" \
                              % z_alt
                        chi2r = chi2r_alt
                    else:
                        mod = self.cb._fit_syst(o_series, z_cand, logN, b, resol, maxfev)
                        #self.systs._add(o_series, z_cand, logN, b, resol)
                        #self.systs._update(cand_mod)
                        msg = "added a %s component at redshift %2.4f" \
                              % (o_series, z_cand)
                        chi2r = chi2r_cand
                    print(prefix, "I'm improving a model at redshift %2.4f "\
                          "(%i/%i): %s (red. chi-squared: %3.2f)...          " \
                          % (o_z, i+1, len(old), msg, chi2r))#, end='\r')
                    chi2r_old = chi2r
                    #bic_old = mod._bic
                    count = 0
                    count_good += 1

                self.cb._update_spec()
                if count_good == 10: break
                if count >= 10: break
                if chi2r<chi2r_thres: break


            if count_good == 0:
                print(prefix, "I've not improved the %s system at redshift "\
                      "%2.4f (%i/%i): I was unable to add useful components."\
                      % (o_series, o_z, i+1, len(old)))
            else:
                print(prefix, "I've improved a model at redshift %2.4f "\
                      "(%i/%i) by adding %i components (red. chi-squared: "\
                      "%3.2f).                                                "\
                      "  " % (o_z, i+1, len(old), count, chi2r))

        #self.systs._clean(chi2r_thres)

        return 0


    def add_syst_slide(self, series='CIV',
                       z_start=0, z_end=6, z_step=2e-4,
                       logN_start=12, logN_end=10, logN_step=-0.2,
                       b_start=8, b_end=9, b_step=1.1,
                       resol=45000, col='y', chi2r_thres=2, maxfev=100):
        """ @brief Fit systems by sliding
        @details Slide a set of Voigt models across a spectrum and fit them
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

        #z_range = np.arange(z_start, z_end, z_step)
        z_range = np.array(self.spec.x/xem_d[series_d[series][0]]-1)
        z_min = np.max([(np.min(self.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_max = np.min([(np.max(self.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
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
            self.systs = SystList(id_start=np.max(systs_old._t['id'])+1)
        else:
            self.systs = SystList()
        chi2a = np.full((len(logN_range),len(b_range),len(z_range)), np.inf)
        #self.corr = np.empty((len(logN_range),len(b_range), 2))
        self.corr = np.empty((len(logN_range)*len(b_range), 4))
        """
        self.corr_logN = logN_range
        self.corr_b = b_range
        self.corr_e = (b_range[0]-b_step*0.5, b_range[-1]+b_step*0.5,
                       logN_range[0]-logN_step*0.5,logN_range[-1]+logN_step*0.5)
        """
        for ilogN, logN in enumerate(logN_range):
            for ib, b in enumerate(b_range):
                icorr = ilogN*len(b_range)+ib
                xm, ym, ym_0, ym_1, ym_2 = self.cb._create_doubl(series, z_mean,
                                                                 logN, b, resol)
                cond_c = 0
                cond_swap_c = 0
                chi2_arr = []
                chi2_0_arr = []
                for iz, z in enumerate(z_range):
                    print(prefix, "I'm testing a %s system (logN=%2.2f, "
                          "b=%2.2f) at redshift %2.4f (%i/%i)..." \
                          % (series, logN, b, z, iz+1, len(z_range)), end='\r')
                    cond, chi2, chi2_0 = \
                        self.cb._test_doubl(xm*(1+z), ym, ym_0, ym_1, ym_2, col)
                    cond_swap, _, _ = \
                        sess.cb._test_doubl(xm*(1+z), ym, ym_0, ym_1, ym_2, col)
                    chi2_arr.append(chi2)
                    chi2_0_arr.append(chi2_0)
                    if cond:
                        chi2a[ilogN, ib, iz] = chi2
                        cond_c += 1
                    if cond_swap:
                        cond_swap_c += 1
                #self.corr[ilogN, ib] = (cond_c, cond_swap_c)

                self.corr[icorr, 0] = logN
                self.corr[icorr, 1] = b
                self.corr[icorr, 2] = cond_c
                self.corr[icorr, 3] = cond_swap_c
                    #1-np.array(cond_swap_c)/np.array(cond_c)
                print(prefix, "I've tested a %s system (logN=%2.2f, "\
                      "b=%2.2f) between redshift %2.4f and %2.4f and found %i "\
                      "coincidences"
                      % (series, logN, b, z_range[0], z_range[-1], cond_c),
                      end='')
                if cond_c != 0:
                    #print(" (estimated correctness=%2.0f%%)."
                    #      % (100*self.corr[ilogN, ib]))
                    print(" (%i in the swapped spectrum)." % cond_swap_c)
                else:
                    print(".")

        # Find candidates choosing the local minima of chi2r at coincidences
        lm = np.logical_and(chi2a < np.inf, detect_local_minima(chi2a))
        chi2m = np.where(lm)

        # ...and then appended
        if hasattr(systs_old, '_t'):
            if len(self.systs._t) > 0:
                self.systs._append(systs_old)
            else:
                self.systs = systs_old

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


        self.cb._update_spec()

        # Save the correctness as a two-entry table - to be modified
        """
        rows = np.array([np.array([logN_range]).T])
        cols = np.array([np.array([np.append([np.nan], b_range)])])
        aprint(rows)
        print(cols)
        print(self.corr)
        self.corr_save = np.concatenate((rows, self.corr), axis=-1)
        self.corr_save = np.concatenate((cols, self.corr_save), axis=0)
        """
        self.corr_save = self.corr
        return 0


    def compl_syst(self, series='CIV', n=100,
                   z_start=0, z_end=6, z_step=1e-2,
                   logN_start=15, logN_end=10, logN_step=-0.2,
                   b_start=8, b_end=9, b_step=1.1,
                   resol=45000, col='y', chi2r_thres=2, maxfev=100):
        """ @brief Estimate completeness
        @details Estimate the completeness of system detection by simulating
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

        z_start, z_end = self.cb._adapt_z(series, z_start, z_end)
        z_range = np.arange(z_start, z_end, z_step)
        logN_range = np.arange(logN_start, logN_end, logN_step)
        b_range = np.arange(b_start, b_end, b_step)

        # Previously fitted systems are left fixed...
        sess = dc(self)

        z_mean = 0.5*(z_start+z_end)
        #self.compl = np.empty((len(z_range), len(logN_range),len(b_range)))
        self.compl = np.empty((len(z_range)*len(logN_range)*len(b_range), 4))
        self.compl_e = (b_range[0]-b_step*0.5, b_range[-1]+b_step*0.5,
                        logN_range[0]-logN_step*0.5,
                        logN_range[-1]+logN_step*0.5)

        z_arr = np.array(self.spec.x/xem_d[series_d[series][0]]-1)
        dz = 2e-4
        compl_sum = 0
        for iz, (zs, ze) in enumerate(zip(z_range[:-1], z_range[1:])):

            for ilogN, logN in enumerate(logN_range):
                for ib, b in enumerate(b_range):
                    icompl = (iz*len(logN_range)+ilogN)*len(b_range)+ib
                    #print(icompl)

                    cond_c = 0
                    n_ok = 0

                    sess.systs = dc(self.systs)
                    sess.cb._append_syst()
                    sess.systs._add(series, z_mean, logN, b+0.5*b_step, resol)

                    xm, ym, ym_0, ym_1, ym_2 = sess.cb._create_doubl(
                        series, z_mean, logN, b, resol)
                    xm_e, ym_e, ym_0_e, ym_1_e, ym_2_e = sess.cb._create_doubl(
                        series, z_mean, logN, b+0.0*b_step, resol)

                    n_fail = 0
                    while n_ok < n:
                        print(prefix, "I'm estimating completeness of %s "
                              "systems (z=[%2.2f, %2.2f], logN=%2.2f, b=%2.2f, "
                              "realization %i/%i)..."
                              % (series, zs, ze, logN, b, n_ok+1, n), end='\r')
                        z_rand = np.random.rand()*(ze-zs)+zs
                        sess.spec = dc(self.spec)
                        fail = sess.cb._apply_doubl(xm_e*(1+z_rand), ym_e)
                        #if not fail or 1==1:
                        n_ok += 1
                        z_round = round(z_rand, 4)
                        z_sel = np.where(np.logical_and(
                            z_arr > z_round-1.5*dz,
                            z_arr < z_round+1.5*dz))
                        for z in z_arr[z_sel]:
                            systs = SystList()
                            cond, chi2, chi2_0 = sess.cb._test_doubl(
                                xm*(1+z), ym, ym_0, ym_1, ym_2, col)
                            if cond and np.abs(z-z_rand)<dz:
                                cond_c += 1
                                break
                        if cond == False:
                            pass
                        #else:
                        #n_fail += 1

                    compl = cond_c/n_ok
                    compl_sum += compl
                    self.compl[icompl, 0] = z_rand
                    self.compl[icompl, 1] = logN
                    self.compl[icompl, 2] = b
                    self.compl[icompl, 3] = compl

        print(prefix, "I've estimated completeness of %s systems "
              "(z=[%2.2f, %2.2f], logN=[%2.2f, %2.2f], b=[%2.2f, %2.2f]); "
              "average was %2.0f%%."
              % (series, z_start, z_end, logN_start, logN_end, b_start, b_end,
                 100*(compl_sum)/np.shape(self.compl)[0]))

        return 0

    def convert_x(self, zem=0, xunit=au.km/au.s):
        """ @brief Convert x axis
        @details Convert the x axis to wavelength or velocity units.
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
        """ @brief Convert y axis
        @details Convert the y axis to flux density units.
        @param yunit Unit of flux density
        @return 0
        """

        yunit = au.Unit(yunit)

        for s in self.seq:
            try:
                getattr(self, s)._convert_y(yunit=yunit)
            except:
                pass
        return 0

    def convolve_gauss(self, std=5, input_col='y', output_col='conv'):
        """@brief Convolve with gaussian
        @details Convolve a spectrum column with a gaussian profile using FFT
        transform.
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
        """ @brief Extract nodes
        @details Extract nodes from a spectrum. Nodes are averages of x and y in
        slices, computed after masking lines.
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
        """ @brief Extract region
        @details Extract a spectral region as a new frame.
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

    def find_peaks(self, col='conv', kind='min', kappa=5.0, append=True):
        """ @brief Find peaks
        @details Find the peaks in a spectrum column. Peaks are the extrema
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

        """
        # Create new session
        kwargs = {'path': self.path, 'name': self.name}
        for s in self.seq:
            try:
                kwargs[s] = getattr(self, s)
            except:
                kwargs[s] = None
            print(kwargs[s])
        if kwargs['spec'] != None:
            new = Session(**kwargs)
        else:
            new = None

        return new
        """
        return 0

    def interp_nodes(self, smooth=0):
        """ @brief Interpolate nodes
        @details Interpolate nodes with a univariate spline to estimate the
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

    def merge_syst(self, series='CIV', v_thres=100):
        """ @brief Merge systems
        @details Merge column densities of systems applying a friend-of-friend
        algorithm.
        @param series Series of transitions
        @param v_thres Velocity threshold for merging (km/s)
        @return 0
        """

        v_thres = float(v_thres) * au.km/au.s

        sel = np.where(self.systs._t['series'] == series)
        self.merge = dc(self.systs._t['z', 'logN', 'dlogN'][sel])
        self.cb._merge_syst(self.merge, v_thres)
        self.merge.sort('z')

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

        # XSHOOTER DAS spectrum
        if instr == 'XSHOOTER' and catg[1:5] == 'SPEC':
            self.spec = format.xshooter_das_spectrum(hdul)

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
                        try:
                            np.savetxt(root+'_compl.dat', self.compl, fmt='%s')
                        except:
                            pass
                        try:
                            np.savetxt(root+'_corr.dat', self.corr, fmt='%s')
                        except:
                            pass
                        try:
                            np.savetxt(root+'_merge.dat', self.merge, fmt='%s')
                        except:
                            pass
                    name = root+'_'+s+'.fits'
                    obj = dc(getattr(self, s))
                    t = obj._t
                    t['x'] = t['x'].to(au.nm)
                    t['xmin'] = t['xmin'].to(au.nm)
                    t['xmax'] = t['xmax'].to(au.nm)
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
        """ @brief Shift from rest frame
        @details Shift x axis from rest frame to the original frame.
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


    def shift_to_rf(self, z=0):
        """ @brief Shift to rest frame
        @details Shift x axis to the rest frame.
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
        """ @brief Simulate a system
        @details Simulate a system by adding Voigt model onto a spectrum.
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
