from .functions import *
from .message import *
from .syst_list import SystList
from .syst_model import SystModel
from .vars import *
from astropy import constants as aconst
from copy import deepcopy as dc
import logging
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import sys

prefix = "[INFO] cookbook_absorbers:"

class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        super(CookbookAbsorbers, self).__init__()
        self._refit_n = 3
        self._chi2rav_thres = 1e-2
        self._chi2r_thres = np.inf
        self._dlogN_thres = np.inf
        self._max_nfev = max_nfev_def


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


    def _mod_ccf(self, mod, ym=None, y=None, verbose=True):
        if eval is None:
            ym = mod.eval(x=mod._xf, params=mod._pars)
        if y is None:
            y = mod._yf

        #w = np.abs(np.gradient(eval))
        #w = w/np.sum(w)
        #ccf = np.correlate(eval, mod._yf)[0]

        #mod._yf = mod.eval(x=mod._xf, params=mod._pars)

        #ccf_same = np.correlate(eval, mod._yf, mode='same')
        #ccf_loc = np.argmax(ccf_same)
        #ccf = np.max(ccf_same)
        ccf = np.dot(ym, y)
        #ccf = np.corrcoef(eval, mod._yf)[1][0]
        #plt.plot(mod._xf, y)
        #plt.plot(mod._xf, ym)
        #plt.scatter(mod._xf[ccf_loc], ccf_same[ccf_loc])
        if verbose:
            logging.info("The data-model CCF is %2.3f." % ccf)
        return ccf


    def _mod_ccf_max(self, mod, vstart=-20, vend=20, dv=1e-2, verbose=True):
        sd = -1*int(np.floor(np.log10(dv)))-1

        xmin = mod._xf[0]
        xmax = mod._xf[-1]
        xmean = 0.5*(xmin+xmax)
        v_shift = np.arange(vstart, vend, dv)
        x_shift = xmean * v_shift/aconst.c.to(au.km/au.s).value
        xstart = xmean * vstart/aconst.c.to(au.km/au.s).value
        xend = xmean * vend/aconst.c.to(au.km/au.s).value
        dx = xmean * dv/aconst.c.to(au.km/au.s).value

        x_osampl = np.arange(xmin+xstart, xmax+xend, dx)
        eval_osampl = mod.eval(x=x_osampl, params=mod._pars)
        eval_ref = mod.eval(x=mod._xf, params=mod._pars)
        #plt.plot(x_osampl,eval_osampl, linewidth=3)
        #x_shift = np.arange(xstart, xend, dx)
        #print(len(x_shift), len(v_shift))
        ccf = []
        grad = np.abs(np.gradient(eval_ref))
        y = (1-mod._yf)*grad/np.sum(grad)
        #plt.plot(mod._xf, mod._yf)
        for xs in x_shift:
            x = x_osampl+xs
            eval = np.interp(mod._xf, x, eval_osampl)
            grad = np.abs(np.gradient(eval))
            ym = (1-eval)#*(grad/np.sum(grad))
            y = (1-mod._yf)#*(grad/np.sum(grad))
            ccf1 = self._mod_ccf(mod, ym, y, verbose=False)
            ccf.append(ccf1)
        amax = np.argmin(ccf)
        #amax = np.argmax(ccf)
        deltax = x_shift[amax]
        deltav = v_shift[amax]
        plt.scatter(xmean+x_shift, ccf/ccf[amax])
        plt.scatter(xmean+x_shift[amax], 1)
        plt.show()
        if verbose:
            logging.info(("I maximized the data model CCF with a shift of "
                          "%."+str(sd)+"e nm (%."+str(sd)+"e km/s)") \
                          % (deltax, deltav))
        return ccf[amax], deltax, deltav


    def _mods_ccf_max(self, vstart, vend, dv):
        systs = self.sess.systs
        for i, m in enum_tqdm(systs._mods_t, len(systs._mods_t),
                              "cookbook_absorbers: Computing CCF"):
            ccf, deltax, deltav = self._mod_ccf_max(m['mod'], vstart, vend, dv,
                                                    verbose=False)
            for i in m['id']:
                w = np.where(systs._t['id']==i)
                systs._t['ccf_deltav'][w] = deltav

        return 0


    def _mods_recreate(self, verbose=True):
        """ Create new system models from a system list """
        spec = self.sess.spec
        systs = self.sess.systs
        #if len(systs._t)==0: return 0
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        #for i,s in enumerate(systs._t):
        for i,s in enum_tqdm(systs._t, len(systs._t),
                             "cookbook_absorbers: Recreating"):
            systs._id = s['id']
            mod = SystModel(spec, systs, z0=s['z0'])
            mod._new_voigt(series=s['series'], z=s['z'], logN=s['logN'],
                           b=s['b'], resol=s['resol'])
            #mod._pars['lines_voigt_%i_z' % i].stderr=s['dz']
            #mod._pars['lines_voigt_%i_logN' % i].stderr=s['dlogN']
            #mod._pars['lines_voigt_%i_b' % i].stderr=s['db']
            #mod._pars['lines_voigt_%i_btur' % i].stderr=0
            #mod._pars['psf_gauss_0_resol'].stderr=0

            self._mods_update(mod)
        mods_n = len(self.sess.systs._mods_t)
        if verbose:
            logging.info("I've recreated %i model%s." \
                         % (mods_n, '' if mods_n==1 else 's'))
        #mod = self.sess.systs._mods_t['mod'][0]
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


    def _resol_update(self, resol):
        mods_t = self.sess.systs._mods_t
        for m in mods_t:
            m['mod']._pars['psf_gauss_0_resol'].value = resol


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


    def _syst_fit(self, mod, verbose=True):
        if self._max_nfev > 0:
            mod._fit(fit_kws={'max_nfev': self._max_nfev})
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
        trans = trans_parse(series)
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


    def _systs_cycle(self, verbose=True):
        chi2rav = np.inf
        chi2rav_old = 0
        chi2r_list, z_list = [], []
        for i,_ in enum_tqdm(range(self._refit_n), self._refit_n,
                              'cookbook_absorbers: Cycling'):
            if chi2rav > self._chi2rav_thres and chi2rav != chi2rav_old:
                if chi2rav < np.inf: chi2rav_old = chi2rav
                chi2r_list, z_list = self._systs_fit(verbose=False)
                if i > 1 and len(chi2r_list)==len(chi2r_list_old):
                    chi2rav = np.mean(np.abs(np.array(chi2r_list)\
                                             -np.array(chi2r_list_old)))
                chi2r_list_old = chi2r_list
                self._systs_reject(verbose=False)
                self._mods_recreate(verbose=False)
            #print(chi2rav, chi2rav_old)
        chi2r_list, z_list = self._systs_fit(verbose=False)
        if verbose and z_list != []:
            logging.info("I've fitted %i model%s." \
                         % (len(self.sess.systs._mods_t), msg_z_range(z_list)))
            if chi2rav < np.inf:
                logging.info("Average chi2r variation after last cycle: %2.4e."\
                             % chi2rav)


    def _systs_fit(self, verbose=True):
        systs = self.sess.systs
        mods_t = systs._mods_t
        if self._max_nfev > 0:
            z_list = []
            chi2r_list = []
            for i,m in enum_tqdm(mods_t, len(mods_t),
                                 "cookbook_absorbers: Fitting"):
            #for i,m in enumerate(mods_t):
                z_list.append(m['z0'])
                self._syst_fit(m['mod'], verbose=False)
                chi2r_list.append(m['mod']._chi2r)

            if verbose:
                logging.info("I've fitted %i model%s." \
                             % (len(mods_t), msg_z_range(z_list)))
        else:
            if verbose:
                logging.info("I've not fitted any model because you choose "
                             "max_nfev=0.")
        return chi2r_list, z_list


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

    """
    def _systs_refit(self, refit_id=[], max_nfev=0):
        systs = self.sess.systs
        mods_t = systs._mods_t
        if max_nfev > 0:
            z_list = []
            mod_list = []
            #chi2r_list = []
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
                #chi2r_list.append(m._chi2r)
            #print(chi2r_list)
            logging.info("I've refitted %i model%s." \
                         % (len(z_list), msg_z_range(z_list)))
            logging.info("I've updated %i system%s in %i model%s." \
                         % (len(systs._t), '' if len(systs._t)==1 else 's',
                            len(mods_t), msg_z_range(z_list)))
        else:
            logging.info("I'm not refitting any system because you choose "
                         "max_nfev=0.")
        return 0
    """

    def _systs_reject(self, verbose=True):
        systs = self.sess.systs
        chi2r_cond = systs._t['chi2r'] > self._chi2r_thres
        """
        relerr_cond = np.logical_or(np.logical_or(
            systs._t['dz'] > dlogN_thres*systs._t['z'],
            systs._t['dlogN'] > dlogN_thres*systs._t['logN']),
            systs._t['db'] > dlogN_thres*systs._t['b'])
        """
        relerr_cond = systs._t['dlogN'] > self._dlogN_thres

        rem = np.where(np.logical_or(chi2r_cond, relerr_cond))[0]
        z_rem = systs._t['z'][rem]
        #refit_id = []
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
            self._systs_remove(rem)#, refit_id)
            """
            for i, r in enum_tqdm(rem, len(rem), "cookbook_absorbers: Removing"):
                t_id = systs._t['id']
                mods_t_id = systs._mods_t['id']
                sel = [t_id[r] in m for m in mods_t_id]
                #if not np.all(np.in1d(mods_t_id[sel][0], t_id[rem])):
                #    refit_id.append(np.setdiff1d(mods_t_id[sel][0], t_id[rem])[0])
            systs._t.remove_rows(rem)
            """
            if verbose:
                logging.info("I've rejected %i mis-identified system%s (%i with a "\
                             "reduced chi2 above %2.2f, %i with relative errors "\
                             "above %2.2f)."
                             % (len(rem), '' if len(rem)==1 else 's',
                                np.sum(chi2r_cond), self._chi2r_thres,
                                np.sum(relerr_cond), self._dlogN_thres))
        return 0 #refit_id


    def _systs_remove(self, rem):#, refit_id):
        systs = self.sess.systs
        for i, r in enum_tqdm(rem, len(rem), "cookbook_absorbers: Removing"):
            t_id = systs._t['id']
            mods_t_id = systs._mods_t['id']
            sel = [t_id[r] in m for m in mods_t_id]
            #if not np.all(np.in1d(mods_t_id[sel][0], t_id[rem])):
            #    refit_id.append(np.setdiff1d(mods_t_id[sel][0], t_id[rem])[0])
        systs._t.remove_rows(rem)


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


    def comp_extract(self, num=1):
        """ @brief Extract systems
        @details Extract systems with less than a given number of components
        @param num Number of components
        @return 0
        """

        try:
            num = int(num)
        except:
            logging.error(msg_param_fail)
            return 0

        out = self.sess.systs
        t_sel = []
        mods_t_sel = []
        for i, m in enumerate(out._mods_t):
            #print(len(m['id']))
            if len(m['id']) > num:
                for id in m['id']:
                    t_sel.append(np.where(out._t['id'] == id)[0])
                mods_t_sel.append(i)
        #print(t_sel)
        out._t.remove_rows(t_sel)
        out._mods_t.remove_rows(mods_t_sel)
        #print(out._t)
        return 0


    def mods_ccf_max(self, vstart=-20, vend=20, dv=1):
        """ @brief Maximize data/model CCF
        @details Slide the system models around their mean wavelength to
        determine the best data/model CCF and the corresponding shift.
        The wavelength range used for sliding is defined in velocity units.
        @param vstart Range start (km/s with respect to mean wavelength)
        @param vend Range end (km/s with respect to mean wavelength)
        @param dv Range step (km/s)
        @return 0
        """

        try:
            vstart = float(vstart)
            vend = float(vend)
            dv = float(dv)
        except:
            logging.error(msg_param_fail)
            return 0

        systs = self.sess.systs
        if 'ccf_deltav' not in systs._t.colnames:
            logging.info("I'm adding column 'ccf_deltav'.")
            systs._t['ccf_deltav'] = np.empty(len(systs._t), dtype=float)

        self._mods_ccf_max(vstart, vend, dv)

        return 0

    def mods_recreate(self):
        """ @brief Recreate the models
        @details Recreate the models from the current system list.
        @return 0
        """

        self._mods_recreate()
        self._spec_update()

        return 0


    def syst_fit(self, num=0, refit_n=0, chi2rav_thres=1e-2,
                 max_nfev=max_nfev_def):
        """ @brief Fit a systems
        @details Fit all Voigt model from a list of systems.
        @param num Number of the system in the list
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            num = int(num)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        mods_t = self.sess.systs._mods_t
        mod = mods_t['mod'][num in mods_t['id']]
        self._syst_fit(mod)
        self._systs_cycle()
        self._spec_update()

        return 0


    def systs_clean(self, chi2r_thres=2.0, dlogN_thres=1.0,
                    max_nfev=max_nfev_def):
        """ @brief Clean system list
        @details Clean systems from a list by rejecting systems with reduced
        chi2 and/or error on column density above a given threshold
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._max_nfev = int(max_nfev)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        self._systs_reject()
        self._mods_recreate()
        self._systs_fit()
        self._spec_update()

        return 0


    def systs_fit(self, refit_n=3, chi2rav_thres=1e-2, max_nfev=max_nfev_def):
        """ @brief Fit systems
        @details Fit all Voigt model from a list of systems.
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        #self._systs_fit()
        self._systs_cycle()
        self._spec_update()

        return 0


### Advanced

    def syst_new(self, series='Ly-a', z=2.0, logN=logN_def, b=b_def,
                 resol=resol_def, chi2r_thres=np.inf, dlogN_thres=np.inf,
                 refit_n=0, chi2rav_thres=1e-2, max_nfev=max_nfev_def):
        """ @brief New system
        @details Add and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @return 0
        """

        try:
            #series = series.replace(';',',')
            z = float(z)
            logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0
        if self._z_off(trans_parse(series), z): return 0

        for s in series.split(';'):
            self._systs_prepare()
        #self._logN_guess(series, z, b, resol)
        #logN = self._syst_guess(series, z)
            mod = self._syst_add(s, z, logN, b, resol)
            if mod is None: return 0
            self._syst_fit(mod)
            self._systs_cycle()
        if self._refit_n == 0:
            self._mods_recreate()
        #refit_id = self._systs_reject(chi2r_thres, dlogN_thres)
        #self._systs_refit(refit_id, max_nfev)
        self._spec_update()

        return 0


    def systs_new_from_lines(self, series='Ly-a', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=np.inf,
                             refit_n=0, chi2rav_thres=1e-2, max_nfev=max_nfev_def,
                             append=True):
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
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @param append Append systems to existing system list
        @return 0
        """

        try:
            #series = series.replace(';',',')
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = float(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0

        for s in series.split(';'):
            z_list, y_list = self._lines_cand_find(s, z_start, z_end, dz)
            z_list, logN_list = self.sess.lines._cand_find2(s, z_start, z_end,
                                                            dz, logN=logN is None)

            if len(z_list) == 0:
                logging.warning("I've found no candidates!")
                return 0

            s_list = [s]*len(z_list)
            resol_list = [resol]*len(z_list)


            self._systs_prepare(append)
            #self._logN_guess(series, z_list[0], b, resol)
            #logN_list = self._systs_guess(series_list, z_list)
            self._systs_add(s_list, z_list, logN_list, resol_list=resol_list)
            #self._systs_fit()
            self._systs_cycle()
            self._spec_update()

        return 0

    def syst_new_from_resids_new(self, series='Ly-a', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=0.5,
                             max_nfev=max_nfev_def, append=True):
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
            series = series.replace(';',',')
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
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
