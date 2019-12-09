from .functions import *
from .message import *
from .syst_list import SystList
from .syst_model import SystModel
from .vars import *
from copy import deepcopy as dc
import logging
import numpy as np
from scipy.interpolate import interp1d
import sys

prefix = "[INFO] cookbook_absorbers:"

class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        pass

    def _lines_cand_find(self, series, z_start, z_end, dz):
        return self.sess.lines._cand_find(series, z_start, z_end, dz)


    def _logN_guess(self, series, z, b, resol):
        spec = dc(self.sess.spec)
        systs = dc(self.sess.systs)
        ynorm_list = []
        logN_list = np.arange(12, 14, 0.1)
        for logN in logN_list:
            mod = SystModel(spec, systs, z0=z)
            mod._new_voigt(series, z, logN, b, resol)
            ynorm_list.append(np.min(mod.eval(x=mod._xs, params=mod._pars)))
        self._guess_f = interp1d(ynorm_list, logN_list-0.5, kind='cubic')


    def _mods_recreate(self, max_nfev=0, verbose=True):
        """ Create new system models from a system list """
        spec = self.sess.spec
        #systs = dc(self.sess.systs)
        systs = self.sess.systs
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        #for i,s in enumerate(systs._t):
        for i,s in enum_tqdm(systs._t, len(systs._t),
                             "cookbook_absorbers: Recreating"):
            systs._id = s['id']
            mod = SystModel(spec, systs, z0=s['z0'])
            mod._new_voigt(series=s['series'], z=s['z0'], logN=s['logN'],
                           b=s['b'], resol=s['resol'])
            self._mods_update(mod)
        #self.sess.systs._mods_t = systs._mods_t
        return 0


    def _mods_update(self, mod, incr=True):
        systs = self.sess.systs
        if mod._group_sel == -1:
            systs._mods_t.add_row([mod._z0, mod, None, []])
        else:
            systs._mods_t[mod._group_sel]['mod'] = mod
        try:
            systs._mods_t[mod._group_sel]['chi2r'] = mod._chi2r
        except:
            systs._mods_t[mod._group_sel]['chi2r'] = np.nan
        systs._mods_t[mod._group_sel]['id'].append(mod._id)

        #systs._id += 1
        systs._id = np.max(systs._t['id'])+1
        return 0


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


    def _syst_add(self, series, z, logN, b, resol, verbose=True):
        systs = self.sess.systs
        spec = self.sess.spec
        if z in systs._t['z0'] \
            and series==systs._t['series'][systs._t['z0']==z]:
            if verbose:
                logging.warning("Redshift %2.4f already exists. Choose another "
                                "one." % z)
            return None

        systs._t.add_row(['voigt_func', series, z, z, None, logN, None, b,
                          None, None, None, systs._id])
        #systs._id = np.max(systs._t['id'])+1
        from .syst_model import SystModel
        mod = SystModel(spec, systs, z0=z)
        mod._new_voigt(series, z, logN, b, resol)

        # When a single system is added, it is stored only on the model table
        self._mods_update(mod, incr=False)
        return mod


    def _syst_fit(self, mod, max_nfev, verbose=True):
        if max_nfev > 0:
            #mod._pars._pretty_print()
            mod._fit(fit_kws={'max_nfev': max_nfev})
            #mod._pars._pretty_print()
            if verbose:
                logging.info("I've fitted 1 model at redshift %2.4f." \
                             % mod._z0)
        else:
            logging.info("I'm not fitting the system because you choose "
                         "max_nfev=0.")

        # When a single system is fitted, it is stored also the system table
        self._systs_update(mod)
        return 0


    def _syst_guess(self, series, z):
        spec = dc(self.sess.spec)
        trans = parse(series)
        x = to_x(z, trans[0]).to(xunit_def).value
        ynorm = np.interp(x, spec.x.to(xunit_def).value,
                          (spec.y/spec._t['cont']).value)
        ynorm = max(0.1, min(0.9, ynorm))
        try:
            ciao
            return max(12, min(14, self._guess_f(ynorm)))
        except:
            #logging.info("I couldn't guess logN for system at redshift %2.4f. "
            #             "I'm using %2.4f instead." % (z, logN_def))
            return logN_def


    def _systs_add(self, series_list, z_list, logN_list=None, b_list=None,
                   resol_list=None, verbose=True):
        if logN_list is None: logN_list = [None]*len(series_list)
        if b_list is None: b_list = [None]*len(series_list)
        if resol_list is None: resol_list = [None]*len(series_list)
        systs_n = 0
        for i, (series, z, logN, b, resol) \
            in enum_tqdm(zip(series_list, z_list, logN_list, b_list, resol_list),
                         len(series_list), "cookbook_absorbers: Adding"):
            #in enumerate(zip(series_list, z_list, logN_list, b_list, resol_list)):
            if logN is None: logN = logN_def
            if b is None: b = b_def
            if resol is None: resol = resol_def
            mod = self._syst_add(series, z, logN, b, resol, False)

            # When many systems are added, they are stored in the system table
            if mod is not None:
                systs_n += 1
                self._systs_update(mod)

        # Improve
        mods_t = self.sess.systs._mods_t
        if verbose:
            logging.info("I've added %i system%s in %i model%s." \
                         % (systs_n, '' if systs_n==1 else 's',
                            len(mods_t), msg_z_range(z_list)))
        return 0


    def _systs_fit(self, max_nfev):
        systs = self.sess.systs
        mods_t = systs._mods_t
        if max_nfev > 0:
            z_list = []
            for i,m in enum_tqdm(mods_t, len(mods_t),
                                 "cookbook_absorbers: Fitting"):
                z_list.append(m['z0'])
                self._syst_fit(m['mod'], max_nfev, verbose=False)

            logging.info("I've fitted %i model%s." \
                         % (len(mods_t), msg_z_range(z_list)))
            self._mods_recreate()#resol)
        else:
            logging.info("I've not fitted any model because you choose "
                         "max_nfev=0.")
        return 0


    def _systs_guess(self, series_list, z_list):
        logN_list = np.array([])
        for series, z in zip(series_list, z_list):
            logN_list = np.append(logN_list, self._syst_guess(series, z))
        return logN_list


    def _systs_prepare(self, append=True):
        systs = self.sess.systs
        if systs != None and len(systs.t) != 0 and append:
            systs._append(SystList(id_start=np.max(systs._t['id'])+1))
        else:
            setattr(self.sess, 'systs', SystList())


    def _systs_refit(self, refit_id=[], max_nfev=0):
        systs = self.sess.systs
        mods_t = systs._mods_t
        if max_nfev > 0:
            z_list = []
            mod_list = []
            #for i,m in enum_tqdm(mods_t, len(mods_t),
            #                     "cookbook_absorbers: Refitting"):
            for i,m in enumerate(mods_t):
                systs_s = [np.where(systs._t['id']==id)[0][0] for id in m['id']]

                # Model has systems in the refit list
                mods_cond = np.any([id in m['id'] for id in refit_id])

                # Systems in the model have different chi2r (resulting from
                # previous fit in separate models)
                chi2r_cond_1 = np.unique(systs._t['chi2r'][systs_s]).size > 1

                # Systems in the model have the same chi2r of another model
                # (resulting from previous fit in the same model)
                chi2r_cond_2 = np.sum([s in systs._t['chi2r'][systs_s] \
                                       for s in systs._t['chi2r']]) > len(systs_s)
                if mods_cond or chi2r_cond_1 or chi2r_cond_2:
                    z_list.append(m['z0'])
                    mod_list.append(m['mod'])
            for i,m in enum_tqdm(mod_list, len(mod_list),
                                 "cookbook_absorbers: Refitting"):
                self._syst_fit(m, max_nfev, verbose=False)
            logging.info("I've refitted %i model%s." \
                         % (len(z_list), msg_z_range(z_list)))
            logging.info("I've updated %i system%s in %i model%s." \
                         % (len(systs._t), '' if len(systs._t)==1 else 's',
                            len(mods_t), msg_z_range(z_list)))
        else:
            logging.info("I'm not refitting any system because you choose "
                         "max_nfev=0.")
        return 0


    def _systs_reject(self, chi2r_thres, dlogN_thres, max_nfev=0):
        systs = self.sess.systs
        chi2r_cond = systs._t['chi2r'] > chi2r_thres
        """
        relerr_cond = np.logical_or(np.logical_or(
            systs._t['dz'] > dlogN_thres*systs._t['z'],
            systs._t['dlogN'] > dlogN_thres*systs._t['logN']),
            systs._t['db'] > dlogN_thres*systs._t['b'])
        """
        relerr_cond = systs._t['dlogN'] > dlogN_thres

        rem = np.where(np.logical_or(chi2r_cond, relerr_cond))[0]
        z_rem = systs._t['z'][rem]
        refit_id = []
        if len(rem) > 0:
            # Check if systems to be rejected are in groups with systems to be
            # preserved, and in case flag the latter for refitting
            """
            for i, r in enum_tqdm(rem, len(rem), "cookbook_absorbers: Rejecting"):
                t_id = systs._t['id']
                mods_t_id = systs._mods_t['id']
                sel = [t_id[r] in m for m in mods_t_id]
                if not np.all(np.in1d(mods_t_id[sel][0], t_id[rem])):
                    refit_id.append(np.setdiff1d(mods_t_id[sel][0], t_id[rem])[0])
            systs._t.remove_rows(rem)
            """
            refit_id = self._syst_remove(self, rem, refit_id)
            logging.info("I've rejected %i mis-identified system%s (%i with a "\
                         "reduced chi2 above %2.2f, %i with relative errors "\
                         "above %2.2f)."
                         % (len(rem), '' if len(rem)==1 else 's',
                            np.sum(chi2r_cond), chi2r_thres,
                            np.sum(relerr_cond), dlogN_thres))
        return refit_id


    def _systs_remove(self, rem, refit_id):
        systs = self.sess.systs
        for i, r in enum_tqdm(rem, len(rem), "cookbook_absorbers: Removing"):
            t_id = systs._t['id']
            mods_t_id = systs._mods_t['id']
            sel = [t_id[r] in m for m in mods_t_id]
            if not np.all(np.in1d(mods_t_id[sel][0], t_id[rem])):
                refit_id.append(np.setdiff1d(mods_t_id[sel][0], t_id[rem])[0])
        systs._t.remove_rows(rem)
        self._mods_recreate()
        return refit_id


    def _systs_update(self, mod, incr=True):
        systs = self.sess.systs
        modw = np.where(mod == systs._mods_t['mod'])[0][0]
        ids = systs._mods_t['id'][modw]
        for i in ids:
            try:
                iw = np.where(systs._t['id']==i)[0][0]
                pref = 'lines_voigt_'+str(i)
                systs._t[iw]['z'] = mod._pars[pref+'_z'].value
                systs._t[iw]['dz'] = mod._pars[pref+'_z'].stderr
                systs._t[iw]['logN'] = mod._pars[pref+'_logN'].value
                systs._t[iw]['dlogN'] = mod._pars[pref+'_logN'].stderr
                systs._t[iw]['b'] = mod._pars[pref+'_b'].value
                systs._t[iw]['db'] = mod._pars[pref+'_b'].stderr
                try:
                    systs._t[iw]['resol'] = mod._pars['psf_gauss_0_resol'].value
                except:
                    systs._t[iw]['resol'] = np.nan
                try:
                    systs._t[iw]['chi2r'] = mod._chi2r
                except:
                    systs._t[iw]['chi2r'] = np.nan
            except:
                pass

        if incr and False:
            systs._id += 1


    def _z_off(self, trans, z):
        for t in trans:
            x = xem_d[t]*(1+z)
            if x < self.sess.spec.x[0] or x > self.sess.spec.x[-1]:
                logging.warning("%s transition at redshift %3.2f is outside "
                                "the spectral range! Please choose a different "
                                "series or redshift." % (t, z))
                return 1
        return 0

### Basic

    def syst_fit(self, num=0, max_nfev=100):
        """ @brief Fit a systems
        @details Fit all Voigt model from a list of systems.
        @param num Number of the system in the list
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            num = int(num)
            max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        self._syst_fit(max_nfev)

        return 0


    def systs_fit(self, max_nfev=100):
        """ @brief Fit systems
        @details Fit all Voigt model from a list of systems.
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        self._systs_fit(max_nfev)
        #self._systs_refit(refit_id, max_nfev)
        self._spec_update()
        return 0


### Advanced

    def syst_new(self, series='Ly-a', z=2.0, logN=logN_def, b=b_def,
                 resol=resol_def, chi2r_thres=np.inf, dlogN_thres=np.inf,
                 max_nfev=100):
        """ @brief New system
        @details Add and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            z = float(z)
            logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            chi2r_thres = float(chi2r_thres)
            dlogN_thres = float(dlogN_thres)
            max_nfev = int(max_nfev)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0
        if self._z_off(parse(series), z): return 0

        self._systs_prepare()
        #self._logN_guess(series, z, b, resol)
        #logN = self._syst_guess(series, z)
        mod = self._syst_add(series, z, logN, b, resol)
        if mod is None: return 0
        self._syst_fit(mod, max_nfev)
        refit_id = self._systs_reject(chi2r_thres, dlogN_thres)
        self._systs_refit(refit_id, max_nfev)
        self._spec_update()

        return 0


    def systs_new_from_lines(self, series='Ly-a', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=np.inf,
                             max_nfev=100, append=True):
        """ @brief New systems from line list
        @details Add and fit Voigt models to a line list, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param logN Guess (logarithmic) column density
        @param b Guess doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
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
            dlogN_thres = float(dlogN_thres)
            resol = None if resol in [None, 'None'] else float(resol)
            max_nfev = int(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0

        z_list, y_list = self._lines_cand_find(series, z_start, z_end, dz)
        z_list, logN_list = self.sess.lines._cand_find2(series, z_start, z_end, dz,
                                                        logN=logN is None)

        if len(z_list) == 0:
            logging.warning("I've found no candidates!")
            return 0

        series_list = [series]*len(z_list)
        resol_list = [resol]*len(z_list)

        self._systs_prepare(append)
        #self._logN_guess(series, z_list[0], b, resol)
        #logN_list = self._systs_guess(series_list, z_list)
        self._systs_add(series_list, z_list, logN_list, resol_list=resol_list)
        self._systs_fit(max_nfev)
        refit_id = self._systs_reject(chi2r_thres, dlogN_thres, max_nfev)
        self._systs_refit(refit_id, max_nfev)
        self._spec_update()

        return 0

    def syst_new_from_resids_new(self, series='Ly-a', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=0.5,
                             max_nfev=100, append=True):
        """ @brief New system from residuals
        @details Add and fit a Voigt model from the strongest residual of
        previously fitted models in the neighborhood.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param logN Guess (logarithmic) column density
        @param b Guess doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
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
            dlogN_thres = float(dlogN_thres)
            resol = float(resol)
            max_nfev = int(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        self.gauss_convolve(std=4, input_col='deabs', output_col='deabs_conv')
        sess = self.peaks_find(col='deabs_conv', kind='min', kappa=3.0, new_sess=True)


        #z_cand =

        return sess
