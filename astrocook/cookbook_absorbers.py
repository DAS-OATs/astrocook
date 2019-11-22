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
