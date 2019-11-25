from .functions import parse
from .vars import *
from copy import deepcopy as dc
import logging
import numpy as np

class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        pass

    def _mods_update(self, resol):
        """ Create new system models from a system list """
        spec = self.sess.spec
        systs = dc(self.sess.systs)
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        for i in range(len(systs._t)):
            s = systs._t[i]
            systs._id = s['id']
            from .syst_model import SystModel
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


    def _syst_add(self, series, z, logN, b, resol):
        systs = self.sess.systs
        spec = self.sess.spec
        systs._t.add_row(['voigt_func', series, z, z, None, logN, None, b, None,
                          None, systs._id])
        from .syst_model import SystModel
        mod = SystModel(spec, systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)
        return mod


    def _systs_reject(self, chi2r_thres, relerr_thres, resol):
        systs = self.sess.systs
        chi2r_cond = systs._t['chi2r'] > chi2r_thres
        relerr_cond = np.logical_or(np.logical_or(
            systs._t['dz'] > relerr_thres*systs._t['z'],
            systs._t['dlogN'] > relerr_thres*systs._t['logN']),
            systs._t['db'] > relerr_thres*systs._t['b'])

        rem = np.where(np.logical_or(chi2r_cond, relerr_cond))[0]
        z_rem = systs._t['z'][rem]
        if len(rem) > 0:
            systs._t.remove_rows(rem)
            logging.info("I've rejected %i mis-identified systems (%i with a "\
                         "reduced chi2 above %2.2f, %i with relative errors "\
                         "above %2.2f)."
                         % (len(rem), np.sum(chi2r_cond), chi2r_thres,
                            np.sum(relerr_cond), relerr_thres))
        self._mods_update(resol)
        return 0


    def syst_new(self, series='Lya', z=2.0, logN=logN_def, b=b_def,
                 resol=resol_def, chi2r_thres=np.inf, relerr_thres=0.1,
                 max_nfev=100):
        """ @brief New system
        @details Add and (optionally) fit a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param logN Guess logarithmic column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param relerr_thres Relative error threshold to accept the fitted model
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)
        chi2r_thres = float(chi2r_thres)
        relerr_thres = float(relerr_thres)
        max_nfev = int(max_nfev)

        if self.sess._z_off(parse(series), z): return 0

        self.sess._systs_prepare()
        systs = self.sess.systs
        mod = self._syst_add(series, z, logN, b, resol)
        if max_nfev > 0:
            mod._fit(fit_kws={'max_nfev': max_nfev})
        else:
            logging.info("I'm not fitting the system because you choose "
                         "max_nfev=0.")
        systs._update(mod)
        self._systs_reject(chi2r_thres, relerr_thres, resol)
        self._spec_update()

        return 0


    def systs_new_from_lines(self, series='Lya', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, relerr_thres=0.1,
                             max_nfev=100):
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
        @param relerr_thres Relative error threshold to accept the fitted model
        @param max_nfev Maximum number of function evaluation
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
        relerr_thres = float(relerr_thres)
        resol = float(resol)
        max_nfev = int(max_nfev)

        z_range, logN_range = self.sess.lines._cand_find(series, z_start, z_end,
                                                         dz, logN=logN is None)

        if len(z_range) == 0:
            logging.warning("I've found no candidates!")
            return 0

        self.sess._systs_prepare()
        systs = self.sess.systs
        for i, (z, l) in enumerate(zip(z_range, logN_range)):
            if l is None: l = logN_def
            mod = self._syst_add(series, z, l, b, resol)
            systs._update(mod)


        mods_t = systs._mods_t
        logging.info("I've added %i %s system(s) in %i model(s) between "
                     "redshift %2.4f and %2.4f." % (len(z_range), series,
                     len(mods_t), z_range[0], z_range[-1]))

        if max_nfev > 0:
            for i,m in enumerate(mods_t):
                print("[INFO] session.add_syst_from_lines: I'm fitting a %s "
                      "model at redshift %2.4f (%i/%i)..."\
                    % (series, m['z0'], i+1, len(mods_t)), end='\r')
                m['mod']._fit(fit_kws={'max_nfev': max_nfev})
                systs._update(m['mod'], mod_t=False)
            try:
                print("[INFO] session.add_syst_from_lines: I've fitted %i %s "
                      "system(s) in %i model(s) between redshift %2.4f and "
                      "%2.4f." \
                      % (len(z_range), series, len(mods_t), z_range[0],
                         z_range[-1]))
            except:
                pass

        self._systs_reject(chi2r_thres, relerr_thres, resol)
        self._spec_update()

        return 0

    def systs_new_from_resids(self, z_start=0, z_end=6, dz=1e-4,
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
                                  peaks.dy, reg._xunit, reg._yunit, reg._meta)
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

                cand_mod = cand.cb._fit_syst(o_series, z_cand, logN, b, resol, maxfev)
                alt_mod = alt.cb._fit_syst('unknown', z_alt, logN, b, resol, maxfev)
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
                        mod = self.sess.cb._fit_syst('unknown', z_alt, logN, b, resol, maxfev)
                        msg = "added an unknown component at wavelength %2.4f" \
                              % z_alt
                        chi2r = chi2r_alt
                    else:
                        mod = self.sess.cb._fit_syst(o_series, z_cand, logN, b, resol, maxfev)
                        msg = "added a %s component at redshift %2.4f" \
                              % (o_series, z_cand)
                        chi2r = chi2r_cand
                    print("[INFO] session.add_syst_from_resids: I'm improving "
                          "a model at redshift %2.4f (%i/%i): %s (red. "
                          "chi-squared: %3.2f)...          " \
                          % (o_z, i+1, len(old), msg, chi2r))#, end='\r')
                    chi2r_old = chi2r
                    count = 0
                    count_good += 1

                self._update_spec()
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
