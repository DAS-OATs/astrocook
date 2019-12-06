from .vars import *
from .functions import *
from .syst_list import SystList
from .syst_model import SystModel
from astropy import units as au
from copy import deepcopy as dc
import logging
import numpy as np

prefix = "[INFO] cookbook_sandbox:"

class CookbookSandbox(object):
    """ Sandbox of utilities to be transferred to other cookbooks
    """

    def __init__(self):
        pass

    def _systs_append(self):
        systs = self.sess.systs
        if systs != None and len(systs.t) != 0:
            #systs._append(SystList(id_start=len(systs._t)))
            systs._append(SystList(id_start=np.max(systs._t['id'])+1))
        else:
            setattr(self.sess, 'systs', SystList())


    def _doubl_apply(self, xm, ym, col='y'):

        spec = self.sess.spec
        y = spec._t[col]
        dy = spec._t['dy']
        s = spec._where_safe

        # If the simulated system is falls by more than a HWHM over a masked
        # line, it is discarded
        ymin = np.min(ym)
        sort = np.argsort(spec.x.to(au.nm))
        ys = np.interp(spec.x.to(au.nm)[sort], xm, ym)
        ysel = np.where(ys < 0.5*(ymin+1))
        #print(0.5*(ymin+1), np.sum(spec._t['lines_mask'][ysel]))
        if np.sum(spec._t['lines_mask'][ysel]) == 0:
            y[s] = ys * y[s]
            dy[s] = np.sqrt(ys) * dy[s]
            #plt.plot(spec.x[s], y[s]/spec.t['cont'][s])
            #plt.plot(xm*10, ym)
            #plt.show()
            return 0
        else:
            #print(spec._t[ysel])
            #plt.plot(xm, ym)
            #plt.plot(spec.x.to(au.nm), spec.y)
            #plt.plot(spec.x.to(au.nm)[spec._t['lines_mask']==1], spec.y[spec._t['lines_mask']==1])
            return 1


    def _doubl_create(self, series='CIV', z_mean=2.0, logN=14, b=10,
                      resol=70000):

        spec = self.sess.spec
        systs = self.sess.systs
        spec._shift_rf(z_mean)
        mod = SystModel(spec, systs, z0=0)
        mod._new_voigt(series, 0, logN, b, resol)
        spec._shift_rf(0.0)
        xm = mod._xm
        hlenm = len(xm)//2
        ym = mod.eval(x=xm, params=mod._pars)
        ym_0 = np.ones(len(xm))
        ym_1 = np.concatenate([ym[:-hlenm], np.ones(hlenm)])
        ym_2 = np.concatenate([np.ones(hlenm), ym[hlenm:]])

        return xm, ym, ym_0, ym_1, ym_2


    def _doubl_test(self, xm, ym, ym_0, ym_1, ym_2, col='y'):
        spec = self.sess.spec
        sort = np.argsort(spec.x.to(au.nm))
        ys = np.interp(xm, spec.x.to(au.nm)[sort], (spec._t[col]/spec._t['cont'])[sort])
        dys = np.interp(xm, spec.x.to(au.nm)[sort], (spec.dy/spec._t['cont'])[sort])
        chi2 = np.sum(((ys-ym)/dys)**2)
        chi2r = chi2/(len(ys)-3)
        chi2_0 = np.sum(((ys-ym_0)/dys)**2)
        chi2_1 = np.sum(((ys-ym_1)/dys)**2)
        chi2_2 = np.sum(((ys-ym_2)/dys)**2)
        fact = 0.7
        if chi2 < fact*np.min([chi2_0-3, chi2_1, chi2_2]):
            return True, chi2, chi2_0
        else:
            return False, chi2, chi2_0

        return 0


    def _mods_update_old(self, resol=70000.0):
        """ Create new system models from a system list """
        spec = self.sess.spec
        systs = dc(self.sess.systs)
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        for i in range(len(systs._t)):
            s = systs._t[i]
            systs._id = s['id']
            mod = SystModel(spec, systs, z0=s['z0'])
            mod._new_voigt(series=s['series'], z=s['z'], logN=s['logN'],
                           b=s['b'], resol=resol)
            systs._update(mod)
        self.sess.systs._mods_t = systs._mods_t


    def _spec_update(self):
        spec = self.sess.spec
        systs = self.sess.systs

        systs._xs = np.array(spec._safe(spec.x).to(au.nm))
        s = spec._where_safe

        y = spec.y
        if 'model' not in spec._t.colnames:
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*y.unit
        if 'deabs' not in spec._t.colnames:
            spec._t['deabs'] = y

        cont = spec._t['cont']
        model = spec._t['model']
        deabs = spec._t['deabs']

        model[s] = cont[s]
        for i, r in enumerate(systs._mods_t):
            mod = r['mod']
            model[s] = mod.eval(x=systs._xs, params=mod._pars) * model[s]
        deabs[s] = cont[s] + y[s] - model[s]


    def _syst_fit_old(self, series='CIV', z=2.0, logN=13.0, b=10.0, resol=70000.0,
                  maxfev=100):

        spec = self.sess.spec
        systs = self.sess.systs
        systs._add(series, z, logN, b, resol)
        mod = SystModel(spec, systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)
        if maxfev > 0:
            mod._fit(fit_kws={'max_nfev': maxfev})
        #plt.plot(mod._xf, mod._yf)
        #plt.show()
        systs._update(mod)
        return mod


    def _syst_merge(self, merge_t, v_thres):

        dv_t = np.array([[ac.c.to(au.km/au.s).value*(z1-z2)/(1+z1)
                          for z1 in merge_t['z']]
                         for z2 in merge_t['z']]) * au.km/au.s
        dv_t[dv_t<=0] = np.inf
        if np.min(dv_t) < v_thres:
            m1, m2 = np.unravel_index(np.argmin(dv_t), np.shape(dv_t))
            z = np.array([merge_t['z'][m1], merge_t['z'][m2]])
            logN = np.array([merge_t['logN'][m1], merge_t['logN'][m2]])
            dlogN = np.array([merge_t['dlogN'][m1], merge_t['dlogN'][m2]])
            z_ave = np.average(z, weights=logN)
            logN_ave = np.log10(np.sum(10**logN))
            dlogN_ave = np.log10(np.sqrt(np.sum(10**(2*dlogN))))
            merge_t.remove_rows([m1, m2])
            merge_t.add_row([z_ave, logN_ave, dlogN_ave])
            self._syst_merge(merge_t, v_thres)

        return 0


    def _syst_mod(self, series='CIV', z=2.0, logN=13.0, b=10.0, resol=70000.0):

        spec = self.sess.spec
        systs = self.sess.systs
        systs._add(series, z, logN, b, resol)
        mod = SystModel(spec, systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)
        systs._update(mod)
        return 0


    def _syst_simul(self, series='Ly_a', z=2.0, logN=13.0, b=10.0,
                    resol=70000.0, col='y'):

        spec = self.sess.spec
        systs = self.sess.systs
        s = spec._where_safe

        systs._add(series, z, logN, b, resol)
        mod = SystModel(spec, systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)

        systs._xs = np.array(spec._safe(spec.x).to(au.nm))
        y = spec._t[col]
        dy = spec._t['dy']
        eval = mod.eval(x=systs._xs, params=mod._pars)

        # If the simulated system is falls by more than a HWHM over a masked
        # line, it is discarded
        ymin = np.min(eval)
        ysel = np.where(eval < 0.5*(ymin+1))
        if np.sum(spec.t['lines_mask'][ysel]) == 0:
            y[s] = eval * y[s]
            dy[s] = np.sqrt(eval) * dy[s]
            return 0
        else:
            return 1


    def _z_adapt(self, series, z_start, z_end):
        spec = self.sess.spec
        z_min = np.max([(np.min(spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_start = max(z_min, z_start)
        z_max = np.min([(np.max(spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_end = min(z_max, z_end)
        return z_start, z_end

    """
    def _fit_mod(self, mod, maxfev=None):
        systs = self.sess.systs
        mod._fit(fit_kws={'max_nfev': maxfev})
        #plt.plot(mod._xf, mod._yf)
        #plt.show()
        #print(systs._t)
        systs._update(mod, mod_t=False)
        #print(systs._t)
    """

    """
    def _rebin(self, dx, xunit):
        spec_in = dc(self.sess.spec)
        spec_in.t.sort('x')
        spec_in._x_convert(xunit=xunit)
        xstart, xend = np.nanmin(spec_in.x), np.nanmax(spec_in.x)
        x = np.arange(xstart.value, xend.value, dx) * xunit
        format = Format()
        xmin, xmax = format._create_xmin_xmax(x)
        im = 0
        iM = 1
        xmin_in = spec_in.xmin[iM].value
        xmax_in = spec_in.xmax[im].value
        y = np.array([])
        dy = np.array([])
        for i, (m, M) \
            in enumerate(tqdm(zip(xmin.value, xmax.value), ncols=100,
                         desc='cookbook._rebin - INFO', total=len(xmin))):
            while xmin_in < M:
                iM += 1
                xmin_in = spec_in.xmin[iM].value
            while xmax_in < m:
                im += 1
                xmax_in = spec_in.xmax[im].value
            ysel = spec_in.y[im:iM+1]
            dysel = spec_in.dy[im:iM+1]
            y = np.append(y, np.average(ysel, weights=1/dysel**2))
            dy = np.append(dy, np.sqrt(np.sum(dysel**2/dysel**4))/np.sum(1/dysel**2))
            #logging.info("I've scanned %i%% of the input spectrum."
            #             % int(100*i/len(xmin)))

        spec_out = Spectrum(x, xmin, xmax, y, dy, xunit=xunit,
                            yunit=spec_in.y.unit, meta=spec_in.meta)
        spec_out._x_convert(xunit=self.sess.spec.x.unit)
        return spec_out
    """

    def syst_simul(self, series='Ly_a', z=2.0, logN=14, b=10, resol=70000,
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

        self._systs_append()
        self._syst_simul(series, z, logN, b, resol, col)

        return 0


    def systs_compl(self, series='CIV', n=100,
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

        z_start, z_end = self._z_adapt(series, z_start, z_end)
        z_range = np.arange(z_start, z_end, z_step)
        logN_range = np.arange(logN_start, logN_end, logN_step)
        b_range = np.arange(b_start, b_end, b_step)

        # Previously fitted systems are left fixed...
        cb = dc(self)
        sess = dc(self.sess)

        z_mean = 0.5*(z_start+z_end)
        #self.sess.compl = np.empty((len(z_range), len(logN_range),len(b_range)))
        self.sess.compl = np.empty((len(z_range)*len(logN_range)*len(b_range), 4))
        self.sess.compl_e = (b_range[0]-b_step*0.5, b_range[-1]+b_step*0.5,
                        logN_range[0]-logN_step*0.5,
                        logN_range[-1]+logN_step*0.5)

        z_arr = np.array(self.sess.spec.x/xem_d[series_d[series][0]]-1)
        dz = 2e-4
        compl_sum = 0
        from .syst_list import SystList
        for iz, (zs, ze) in enumerate(zip(z_range[:-1], z_range[1:])):

            for ilogN, logN in enumerate(logN_range):
                for ib, b in enumerate(b_range):
                    icompl = (iz*len(logN_range)+ilogN)*len(b_range)+ib
                    #print(icompl)

                    cond_c = 0
                    n_ok = 0

                    sess.systs = dc(self.sess.systs)
                    self._systs_append()
                    sess.systs._add(series, z_mean, logN, b+0.5*b_step, resol)

                    xm, ym, ym_0, ym_1, ym_2 = cb._doubl_create(
                        series, z_mean, logN, b, resol)
                    xm_e, ym_e, ym_0_e, ym_1_e, ym_2_e = cb._doubl_create(
                        series, z_mean, logN, b+0.0*b_step, resol)

                    n_fail = 0
                    while n_ok < n:
                        print(prefix, "I'm estimating completeness of %s "
                              "systems (z=[%2.2f, %2.2f], logN=%2.2f, b=%2.2f, "
                              "realization %i/%i)..."
                              % (series, zs, ze, logN, b, n_ok+1, n), end='\r')
                        z_rand = np.random.rand()*(ze-zs)+zs
                        sess.spec = dc(self.sess.spec)
                        fail = self._doubl_apply(xm_e*(1+z_rand), ym_e)
                        #if not fail or 1==1:
                        n_ok += 1
                        z_round = round(z_rand, 4)
                        z_sel = np.where(np.logical_and(
                            z_arr > z_round-1.5*dz,
                            z_arr < z_round+1.5*dz))
                        for z in z_arr[z_sel]:
                            systs = SystList()
                            cond, chi2, chi2_0 = cb._doubl_test(
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
                    self.sess.compl[icompl, 0] = z_rand
                    self.sess.compl[icompl, 1] = logN
                    self.sess.compl[icompl, 2] = b
                    self.sess.compl[icompl, 3] = compl

        logging.info("I've estimated completeness of %s systems "
              "(z=[%2.2f, %2.2f], logN=[%2.2f, %2.2f], b=[%2.2f, %2.2f]); "
              "average was %2.0f%%."
              % (series, z_start, z_end, logN_start, logN_end, b_start, b_end,
                 100*(compl_sum)/np.shape(self.sess.compl)[0]))

        return 0


    def systs_merge(self, series='CIV', v_thres=100):
        """ @brief Merge systems
        @details Merge column densities of systems applying a friend-of-friend
        algorithm.
        @param series Series of transitions
        @param v_thres Velocity threshold for merging (km/s)
        @return 0
        """

        v_thres = float(v_thres) * au.km/au.s

        sel = np.where(self.sess.systs._t['series'] == series)
        self.sess.merge = dc(self.sess.systs._t['z', 'logN', 'dlogN'][sel])
        self._syst_merge(self.sess.merge, v_thres)
        self.sess.merge.sort('z')

        return 0


    def systs_new_from_resids(self, z_start=0, z_end=6, dz=1e-4,
                              resol=resol_def, logN=11, b=5, chi2r_thres=1.0,
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
        resol = None if resol in [None, 'None'] else float(resol)#resol = float(resol)
        logN = float(logN)
        b = float(b)
        chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        systs = self.sess.systs

        # Select systems by redshift range and reduced chi2 threshold
        cond_z =  np.logical_and(systs._mods_t['z0'] > z_start,
                                 systs._mods_t['z0'] < z_end)
        cond_chi2r = np.logical_or(systs._mods_t['chi2r'] > chi2r_thres,
                                   np.isnan(systs._mods_t['chi2r']))
        old = systs._mods_t[np.where(np.logical_and(cond_z, cond_chi2r))]

        for i, o in enumerate(old):
            o_id = o['id'][0]
            o_series = systs._t[systs._t['id'] == o_id]['series'][0]
            o_z = np.array(systs._t['z'][systs._t['id']==o_id])[0]

            chi2r_old = np.inf
            count = 0
            count_good = 0

            while True:

                spec = dc(self.sess.spec)
                spec._gauss_convolve(std=2, input_col='deabs', verb=False)
                try:
                    reg_x = systs._mods_t['mod'][i]._xm
                except:
                    break

                reg_xmin = np.interp(reg_x, spec.x.to(au.nm), spec.xmin.to(au.nm))
                reg_xmax = np.interp(reg_x, spec.x.to(au.nm), spec.xmax.to(au.nm))
                reg_y = np.interp(reg_x, spec.x.to(au.nm), spec.t['conv']-spec.t['cont'])
                reg_dy = np.interp(reg_x, spec.x.to(au.nm), spec.dy)
                from .spectrum import Spectrum
                reg = Spectrum(reg_x, reg_xmin, reg_xmax, reg_y, reg_dy)
                peaks = reg._peaks_find(col='y')#, mode='wrap')

                from .line_list import LineList
                resids = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.y,
                                  peaks.dy, xunit=reg._xunit, yunit=reg._yunit, meta=reg._meta)
                z_cand = resids._syst_cand(o_series, z_start, z_end, dz,
                                           single=True)
                z_alt = resids._syst_cand('unknown', 0, np.inf, dz, single=True)

                # If no residuals are found, add a system at the init. redshift
                if z_cand == None:
                    z_cand = o_z

                if z_alt == None:
                    z_alt = (1.+o_z)\
                            *xem_d[series_d[o_series][0]].to(au.nm).value

                if count == 0:
                    t_old, mods_t_old = self.sess.systs._freeze()
                from .syst_list import SystList
                systs._append(SystList(id_start=np.max(self.sess.systs._t['id'])+1),
                              unique=False)
                cand = dc(self.sess)
                alt = dc(self.sess)
                cand_cb = dc(self)
                alt_cb = dc(self)

                cand_mod = cand_cb._syst_fit_old(o_series, z_cand, logN, b, resol, maxfev)
                alt_mod = alt_cb._syst_fit_old('unknown', z_alt, logN, b, resol, maxfev)
                self.sess.systs._mods_t = dc(mods_t_old)
                chi2r_cand = cand_mod._chi2r
                chi2r_alt = alt_mod._chi2r
                if cand_mod._chi2r>=chi2r_old and alt_mod._chi2r>= chi2r_old:
                    self.sess.systs._unfreeze(t_old, mods_t_old)
                    count += 1
                    break
                else:
                    t_old, mods_t_old = self.sess.systs._freeze()
                    if chi2r_cand > chi2r_alt:#*2:#1.1:# and count > 3:
                        mod = self._syst_fit_old('unknown', z_alt, logN, b, resol, maxfev)
                        msg = "added an unknown component at wavelength %2.4f" \
                              % z_alt
                        chi2r = chi2r_alt
                    else:
                        mod = self._syst_fit_old(o_series, z_cand, logN, b, resol, maxfev)
                        msg = "added a %s component at redshift %2.4f" \
                              % (o_series, z_cand)
                        chi2r = chi2r_cand
                    print(prefix, "I'm improving "
                          "a model at redshift %2.4f (%i/%i): %s (red. "
                          "chi-squared: %3.2f)...          " \
                          % (o_z, i+1, len(old), msg, chi2r), end='\r')
                    chi2r_old = chi2r
                    count = 0
                    count_good += 1

                self._spec_update()
                if count_good == 10: break
                if count >= 10: break
                if chi2r<chi2r_thres: break


            if count_good == 0:
                logging.warning("I've not improved the %s system at redshift "\
                                "%2.4f (%i/%i): I was unable to add useful "\
                                "components."\
                                % (o_series, o_z, i+1, len(old)))
            else:
                logging.info("I've improved a model at redshift %2.4f (%i/%i) "
                             "by adding %i components (red. chi-squared: "\
                             "%3.2f).                                                "\
                             "  " % (o_z, i+1, len(old), count, chi2r))

        return 0

    def systs_new_from_slide(self, series='CIV',
                        z_start=0, z_end=6, z_step=2e-4,
                        logN_start=12, logN_end=10, logN_step=-0.2,
                        b_start=8, b_end=9, b_step=1.1,
                        resol=resol_def, col='y', chi2r_thres=2, maxfev=100):
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
        z_range = np.array(self.sess.spec.x/xem_d[series_d[series][0]]-1)
        z_min = np.max([(np.min(self.sess.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_max = np.min([(np.max(self.sess.spec.x.to(au.nm))/xem_d[t]).value-1.0 \
                        for t in series_d[series]])
        z_range = z_range[np.where(np.logical_and(z_range > z_min,
                                                  z_range < z_max))]
        z_mean = 0.5*(z_min+z_max)
        logN_range = np.arange(logN_start, logN_end, logN_step)
        b_range = np.arange(b_start, b_end, b_step)

        # Create x-swapped spectrum to monitor correctness
        sess = dc(self.sess)
        cb = dc(self)
        sess.spec._t['x'] = sess.spec._t['x'][::-1]

        # Previously fitted systems are left fixed...
        systs_old = dc(self.sess.systs)

        from .syst_list import SystList
        if hasattr(systs_old, '_t'):
            self.sess.systs = SystList(id_start=np.max(systs_old._t['id'])+1)
        else:
            self.sess.systs = SystList()
        chi2a = np.full((len(logN_range),len(b_range),len(z_range)), np.inf)
        #self.sess.corr = np.empty((len(logN_range),len(b_range), 2))
        self.sess.corr = np.empty((len(logN_range)*len(b_range), 4))
        """
        self.sess.corr_logN = logN_range
        self.sess.corr_b = b_range
        self.sess.corr_e = (b_range[0]-b_step*0.5, b_range[-1]+b_step*0.5,
                       logN_range[0]-logN_step*0.5,logN_range[-1]+logN_step*0.5)
        """
        for ilogN, logN in enumerate(logN_range):
            for ib, b in enumerate(b_range):
                icorr = ilogN*len(b_range)+ib
                xm, ym, ym_0, ym_1, ym_2 = self._doubl_create(series, z_mean,
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
                        self._doubl_test(xm*(1+z), ym, ym_0, ym_1, ym_2, col)
                    cond_swap, _, _ = \
                        cb._doubl_test(xm*(1+z), ym, ym_0, ym_1, ym_2, col)
                    chi2_arr.append(chi2)
                    chi2_0_arr.append(chi2_0)
                    if cond:
                        chi2a[ilogN, ib, iz] = chi2
                        cond_c += 1
                    if cond_swap:
                        cond_swap_c += 1
                #self.sess.corr[ilogN, ib] = (cond_c, cond_swap_c)

                self.sess.corr[icorr, 0] = logN
                self.sess.corr[icorr, 1] = b
                self.sess.corr[icorr, 2] = cond_c
                self.sess.corr[icorr, 3] = cond_swap_c
                    #1-np.array(cond_swap_c)/np.array(cond_c)
                print(prefix, "I've tested a %s system (logN=%2.2f, "\
                      "b=%2.2f) between redshift %2.4f and %2.4f and found %i "\
                      "coincidences"
                      % (series, logN, b, z_range[0], z_range[-1], cond_c),
                      end='')
                if cond_c != 0:
                    #print(" (estimated correctness=%2.0f%%)."
                    #      % (100*self.sess.corr[ilogN, ib]))
                    print(" (%i in the swapped spectrum)." % cond_swap_c)
                else:
                    print(".")

        # Find candidates choosing the local minima of chi2r at coincidences
        lm = np.logical_and(chi2a < np.inf, detect_local_minima(chi2a))
        chi2m = np.where(lm)

        # ...and then appended
        if hasattr(systs_old, '_t'):
            if len(self.sess.systs._t) > 0:
                self.sess.systs._append(systs_old)
            else:
                self.sess.systs = systs_old

        logging.info("I've selected %i candidates among the coincidences."\
                     % len(chi2m[0]))
        if maxfev > 0:
            for i in range(len(chi2m[0])):
                z = z_range[chi2m[2][i]]
                logN = logN_range[chi2m[0][i]]
                b = b_range[chi2m[1][i]]
                self._syst_fit_old(series, z, logN, b, resol, maxfev)
                print(prefix, "I've fitted a %s system at redshift %2.4f "
                      "(%i/%i)..." % (series, z, i+1, len(chi2m[0])), end='\r')
            if len(chi2m[0]) > 0:
                print(prefix, "I've fitted %i %s systems between redshift "
                      "%2.4f and %2.4f."
                      % (len(chi2m[0]), series, z_range[chi2m[2][0]],
                         z_range[chi2m[2][-1]]))


        self._spec_update()

        # Save the correctness as a two-entry table - to be modified
        """
        rows = np.array([np.array([logN_range]).T])
        cols = np.array([np.array([np.append([np.nan], b_range)])])
        aprint(rows)
        print(cols)
        print(self.sess.corr)
        self.sess.corr_save = np.concatenate((rows, self.sess.corr), axis=-1)
        self.sess.corr_save = np.concatenate((cols, self.sess.corr_save), axis=0)
        """
        self.sess.corr_save = self.sess.corr
        return 0

    """
    self civ_full(self):

        sess_start = self.sess
        if sess_start.spec.meta['object'] == 'J2123-0050':
            sess = sess_start.region_extract(xmin=xmin, xmax=xmax)
        else:
            sess = sess_start
        sess.gauss_convolve(std=10)
        sess.peaks_find(kappa=3.0)
        sess.lines._t.remove_rows(sess.lines.y == 0)
        if np.mean(sess.spec._t['y'])<1 and np.std(sess.spec._t['y'])<1:
            sess.spec._t['cont'] = [1] * len(sess.spec._t)*sess.spec.y.unit
        if 'cont' not in sess.spec._t.colnames:
            sess.nodes_extract(delta_x=1000)
            sess.nodes_interp()
        sess_reg = sess.region_extract(xmin=xmin, xmax=xmax)
        self._gui._panel_sess._on_add(sess_reg, open=False)
        sess_center = dc(sess_reg)
        sess_center.add_syst_from_lines(z_end=20, maxfev=10)#series='unknown')
        sess_reg.lines.t['x'] = (1+sess_center.systs.t['z'])\
                                *xem_d['Ly_a'].to(sess_reg.spec.x.unit)
        sess_reg.lines.t['logN'] = sess_center.systs.t['logN']
        #sess_reg.add_syst_from_lines(series='SiII', logN=None, b=20.0,
        #                             dz=5e-5, z_end=zem, maxfev=10)
        #sess_reg.add_syst_from_lines(series='SiIV', logN=None, b=20.0,
        #                             dz=5e-5, z_end=zem, maxfev=10)
        sess_reg.add_syst_from_lines(series='CIV', logN=None, b=20.0,
                                     dz=5e-5, z_end=zem, maxfev=10)
        #sess_reg.add_syst_from_lines(series='FeII', logN=None, b=20.0,
        #                             dz=5e-5, z_end=zem, maxfev=10)
        #sess_reg.add_syst_from_lines(series='MgII', logN=None, b=20.0,
        #                             dz=5e-5, z_end=zem, maxfev=10)
        sess_reg.add_syst_from_resids(chi2r_thres=2.0, logN=13.0, b=10.0,
                                      maxfev=10)
        sess_reg.compl_syst(n=10)#, z_start=2.128, z_end=2.1372)
        sess_reg.add_syst_slide(col='deabs')#, z_start=1.6, z_end=1.61)
        sess_reg.syst_merge()
        self._gui._refresh()
        sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_'
                      +datetime.date.today().isoformat()+'.xxx')
        sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_latest.xxx')
        time_end = datetime.datetime.now()
        print("%s; computation time: %s" \
              % (datetime.datetime.now(), time_end-time_start))
    """
