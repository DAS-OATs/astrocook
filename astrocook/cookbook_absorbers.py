from .functions import parse
from .message import *
from .syst_list import SystList
from .vars import *
from copy import deepcopy as dc
import logging
import numpy as np
import sys
from tqdm import tqdm

prefix = "[INFO] cookbook_absorbers:"

class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        pass

    def _mods_recreate(self, resol, max_nfev=0, verbose=True):
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
            #systs._update(mod, t=False)
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
        if incr:
            systs._id += 1



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

        # When a single system is added, it is stored only on the model table
        self._mods_update(mod, incr=False)
        return mod


    def _syst_fit(self, mod, max_nfev):
        if max_nfev > 0:
            mod._fit(fit_kws={'max_nfev': max_nfev})
        else:
            logging.info("I'm not fitting the system because you choose "
                         "max_nfev=0.")

        # When a single system is fitted, it is stored also the system table
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

            # When many systems are added, they are stored in the system table
            self._systs_update(mod)

        # Improve
        mods_t = self.sess.systs._mods_t
        if verbose:
            logging.info("I've added %i system%s in %i model%s." \
                         % (len(z_list), '' if len(z_list)==1 else 's',
                            len(mods_t), msg_z_range(z_list)))
        return 0


    def _systs_fit(self, resol, max_nfev):
        systs = self.sess.systs
        mods_t = systs._mods_t
        if max_nfev > 0:
            z_list = []
            for i,m in enumerate(
                tqdm(mods_t, ncols=120, total=len(mods_t),
                     desc="[INFO] cookbook_absorbers: Fitting")):
                z_list.append(m['z0'])
                self._syst_fit(m['mod'], max_nfev)

                # When many systems are fitted, they are stored again in the
                # system table after fitting
                self._systs_update(m['mod'])
            logging.info("I've fitted %i model%s." \
                         % (len(mods_t), msg_z_range(z_list)))
            self._mods_recreate(resol)
        else:
            logging.info("I've not fitted any model because you choose "
                         "max_nfev=0.")
        return 0


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
            for i,m in enumerate(
                tqdm(mods_t, ncols=120, total=len(mods_t),
                     desc="[INFO] cookbook_absorbers: Refitting")):
                systs_s = [np.where(systs._t['id']==id)[0][0] for id in m['id']]
                mods_s = np.any([id in m['id'] for id in refit_id])
                if np.unique(systs._t['chi2r'][systs_s]).size > 1 or mods_s:
                    z_list.append(m['z0'])
                    self._syst_fit(m['mod'], max_nfev)
                    self._systs_update(m['mod'])
            logging.info("I've refitted %i model%s." \
                         % (len(z_list), msg_z_range(z_list)))
        else:
            logging.info("I'm not refitting any system because you choose "
                         "max_nfev=0.")
        #self._mods_recreate(resol)
        return 0


    def _systs_reject(self, chi2r_thres, dlogN_thres, resol, max_nfev=0):
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
                            np.sum(relerr_cond), dlogN_thres))
        #self._mods_recreate(resol, refit_id, max_nfev)
        self._mods_recreate(resol)
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
                    systs._t[iw]['chi2r'] = mod._chi2r
                except:
                    systs._t[iw]['chi2r'] = np.nan
            except:
                pass

        if incr:
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

################################################################################

    def syst_new(self, series='Lya', z=2.0, logN=logN_def, b=b_def,
                 resol=resol_def, chi2r_thres=np.inf, dlogN_thres=0.5,
                 max_nfev=100):
        """ @brief New system
        @details Add and (optionally) fit a Voigt model for a system.
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
            resol = float(resol)
            chi2r_thres = float(chi2r_thres)
            dlogN_thres = float(dlogN_thres)
            max_nfev = int(max_nfev)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        if self._z_off(parse(series), z): return 0

        self._systs_prepare()
        mod = self._syst_add(series, z, logN, b, resol)
        self._syst_fit(mod, max_nfev)
        refit_id = self._systs_reject(chi2r_thres, dlogN_thres, resol)
        self._systs_refit(refit_id, max_nfev)
        self._spec_update()

        return 0


    def systs_new_from_lines(self, series='Lya', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=0.5,
                             max_nfev=100, append=True):
        """ @brief Fit systems from line list
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
        self._systs_add([series]*len(z_list), z_list, logN_list)
        self._systs_fit(resol, max_nfev)
        refit_id = self._systs_reject(chi2r_thres, dlogN_thres, resol, max_nfev)
        self._systs_refit(refit_id, max_nfev)
        self._spec_update()

        return 0
