from .functions import parse
from .message import *
from .syst_list import SystList
from .vars import *
from copy import deepcopy as dc
import logging
import numpy as np
import sys

prefix = "[INFO] cookbook_absorbers:"

class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        pass

    def _mods_update(self, resol, refit_id=[], max_nfev=0, verbose=True):
        """ Create new system models from a system list """
        spec = self.sess.spec
        #systs = dc(self.sess.systs)
        systs = self.sess.systs
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        for i in range(len(systs._t)):
            s = systs._t[i]
            systs._id = s['id']
            from .syst_model import SystModel
            mod = SystModel(spec, systs, z0=s['z0'])
            mod._new_voigt(series=s['series'], z=s['z0'], logN=s['logN'],
                           b=s['b'], resol=resol)
            if systs._id in refit_id and max_nfev>0:
                if verbose:
                    print(prefix, "I'm refitting a model at redshift %2.4f, "
                          "starting from the parameters of the previous fit..."
                          % s['z'], end='\r')
                mod._fit(fit_kws={'max_nfev': max_nfev})
                systs._update(mod)
            else:
                systs._update(mod, t=False)
        #self.sess.systs._mods_t = systs._mods_t
        return 0


    def _lines_cand_find(self, series, z_start, z_end, dz, logN):
        return self.sess.lines._cand_find(series, z_start, z_end, dz, logN)


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
        return 0


    def _syst_add(self, series, z, logN, b, resol):
        systs = self.sess.systs
        spec = self.sess.spec
        systs._t.add_row(['voigt_func', series, z, z, None, logN, None, b, None,
                          None, systs._id])
        from .syst_model import SystModel
        mod = SystModel(spec, systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)
        return mod


    def _syst_fit(self, mod, max_nfev):
        if max_nfev > 0:
            mod._fit(fit_kws={'max_nfev': max_nfev})
        else:
            logging.info("I'm not fitting the system because you choose "
                         "max_nfev=0.")
        self._systs_update(mod)
        return 0


    def _systs_add(self, series_list, z_list, logN_list=None, b_list=None,
                   resol_list=None, verbose=True):
        if logN_list is None: logN_list = [None]*len(series_list)
        if b_list is None: b_list = [None]*len(series_list)
        if resol_list is None: resol_list = [None]*len(series_list)
        for i, (series, z, logN, b, resol) \
            in enumerate(zip(series_list, z_list, logN_list, b_list, resol_list)):
            if logN is None: logN = logN_def
            if b is None: b = b_def
            if resol is None: resol = resol_def
            mod = self._syst_add(series, z, logN, b, resol)
            self._systs_update(mod)

        # Improve
        mods_t = self.sess.systs._mods_t
        if verbose:
            logging.info("I've added %i system%s in %i model%s between "
                         "redshift %2.4f and %2.4f." \
                         % (len(z_list), '' if len(z_list)==1 else 's',
                            len(mods_t), '' if len(mods_t)==1 else 's',
                            z_list[0], z_list[-1]))
        return 0


    def _systs_fit(self, max_nfev, verbose=True):
        if max_nfev > 0:
            systs = self.sess.systs
            mods_t = systs._mods_t
            z_list = []
            for i,m in enumerate(mods_t):
                z_list.append(m['z0'])
                if verbose:
                    print(prefix, "I'm fitting a model at redshift %2.4f "
                          "(%i/%i)..."
                          % (m['z0'], i+1, len(mods_t)), end='\r')
                m['mod']._fit(fit_kws={'max_nfev': max_nfev})
                systs._update(m['mod'], mod_t=False)
            if verbose:
                print(prefix, "I've fitted %i model%s between redshift %2.4f "
                      "and %2.4f." \
                      % (len(mods_t), '' if len(mods_t)==1 else 's',
                         np.min(z_list), np.max(z_list)))
        else:
            logging.info("I'm not fitting any system because you choose "
                         "max_nfev=0.")
        return 0


    def _systs_prepare(self, append=True):
        systs = self.sess.systs
        if systs != None and len(systs.t) != 0 and append:
            systs._append(SystList(id_start=np.max(systs._t['id'])+1))
        else:
            setattr(self.sess, 'systs', SystList())


    def _systs_refit(self, refit_id=[], max_nfev=0, verbose=True):
        systs = self.sess.systs
        mods_t = systs._mods_t
        if max_nfev > 0:
            z_refit = []
            for m in mods_t:
                systs_s = [np.where(systs._t['id']==id)[0][0] for id in m['id']]
                mods_s = np.any([id in m['id'] for id in refit_id])
                if np.unique(systs._t['chi2r'][systs_s]).size > 1 or mods_s:
                    if verbose:
                        print(prefix, "I'm refitting a model at redshift "
                              "%2.4f, starting from the result of the previous "
                              "fit..." % m['z0'], end='\r')
                    z_refit.append(m['z0'])
                    m['mod']._fit(fit_kws={'max_nfev': max_nfev})
                    systs._update(m['mod'])
            if verbose:
                print(prefix, "I've refitted %i model%s between redshift %2.4f "
                      "and %2.4f." \
                      % (len(z_refit), '' if len(z_refit)==1 else 's',
                         np.min(z_refit), np.max(z_refit)))
        else:
            logging.info("I'm not refitting any system because you choose "
                         "max_nfev=0.")
        return 0


    def _systs_reject(self, chi2r_thres, relerr_thres, resol, max_nfev=0):
        systs = self.sess.systs
        chi2r_cond = systs._t['chi2r'] > chi2r_thres
        relerr_cond = np.logical_or(np.logical_or(
            systs._t['dz'] > relerr_thres*systs._t['z'],
            systs._t['dlogN'] > relerr_thres*systs._t['logN']),
            systs._t['db'] > relerr_thres*systs._t['b'])

        rem = np.where(np.logical_or(chi2r_cond, relerr_cond))[0]
        z_rem = systs._t['z'][rem]
        refit_id = []
        if len(rem) > 0:
            # Check if systems to be rejected are in groups with systems to be
            # preserved, and in case flag the latter for refitting
            for r in rem:
                t_id = systs._t['id']
                mods_t_id = systs._mods_t['id']
                sel = [t_id[r] in m for m in mods_t_id]
                if not np.all(np.in1d(mods_t_id[sel][0], t_id[rem])):
                    refit_id.append(np.setdiff1d(mods_t_id[sel][0], t_id[rem])[0])
            systs._t.remove_rows(rem)
            logging.info("I've rejected %i mis-identified system%s (%i with a "\
                         "reduced chi2 above %2.2f, %i with relative errors "\
                         "above %2.2f)."
                         % (len(rem), '' if len(rem)==1 else 's',
                            np.sum(chi2r_cond), chi2r_thres,
                            np.sum(relerr_cond), relerr_thres))
        #self._mods_update(resol, refit_id, max_nfev)
        return refit_id


    def _systs_update(self, mod):
        self.sess.systs._update(mod)


    def _z_off(self, trans, z):
        for t in trans:
            x = xem_d[t]*(1+z)
            if x < self.sess.spec.x[0] or x > self.sess.spec.x[-1]:
                logging.warning("%s transition at redshift %3.2f is outside "
                                "the spectral range! Please choose a different "
                                "series or redshift." % (t, z))
                return 1
        return 0

################################################################################

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

        try:
            z = float(z)
            logN = float(logN)
            b = float(b)
            resol = float(resol)
            chi2r_thres = float(chi2r_thres)
            relerr_thres = float(relerr_thres)
            max_nfev = int(max_nfev)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        if self._z_off(parse(series), z): return 0

        self._systs_prepare()
        mod = self._syst_add(series, z, logN, b, resol)
        self._syst_fit(mod, max_nfev)
        self._systs_reject(chi2r_thres, relerr_thres, resol)
        self._spec_update()

        return 0


    def systs_new_from_lines(self, series='Lya', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, relerr_thres=0.1,
                             max_nfev=100, append=True):
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
        @param append Append systems to existing system list
        @return 0
        """

        try:
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
            append = append == 'True'
        except ValueError:
            logging.error(msg_param_fail)
            return 0


        z_list, logN_list = self._lines_cand_find(series, z_start, z_end, dz,
                                                   logN=logN is None)

        if len(z_list) == 0:
            logging.warning("I've found no candidates!")
            return 0

        self._systs_prepare(append)
        #"""
        self._systs_add([series]*len(z_list), z_list, logN_list)
        self._systs_fit(max_nfev)
        """
        for zi, logNi, in zip(z_list, logN_list):
            print(zi)
            mod = self._syst_add(series, zi, logNi, b, resol)
            self._syst_fit(mod, max_nfev)
        """
        #refit_id = self._systs_reject(chi2r_thres, relerr_thres, resol, max_nfev)
        #print(refit_id)
        #refit_id = []
        self._mods_update(resol)
        #self._systs_refit(refit_id, max_nfev)
        self._spec_update()

        return 0
