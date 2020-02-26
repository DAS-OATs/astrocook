from .functions import *
from .message import *
from .syst_list import SystList
from .syst_model import SystModel
from .vars import *
from astropy import constants as aconst
from astropy import table as at
from astropy import units as au
from copy import deepcopy as dc
import logging
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import sys

prefix = "[INFO] cookbook_absorbers:"

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


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


    def _mod_ccf(self, mod, ym=None, y=None, verbose=True, plot=False):
        if ym is None:
            ym = mod.eval(x=mod._xf, params=mod._pars)
        if y is None:
            y = mod._yf

        ccf = np.dot(ym, y)
        if plot:
            plt.plot(mod._xf, ym)
        if verbose:
            logging.info("The data-model CCF is %2.3f." % ccf)
        return ccf


    def _mod_ccf_max(self, mod, vstart=-5, vend=5, dv=1e-2, weight=True,
                     verbose=True):
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
        eval_osampl = 1-mod.eval(x=x_osampl, params=mod._pars)
        ccf = []

        y = (1-mod._yf)
        if weight:
            #w = np.abs(np.gradient(eval_osampl))
            #eval_osampl = eval_osampl * w/np.sum(w)*len(w)
            y = y*mod._wf

        #y = (1-mod._yf)#*grad/np.sum(grad)

        for i, xs in enumerate(x_shift):
            plot = False
            x = x_osampl+xs
            digitized = np.digitize(x, mod._xf)
            ym = [eval_osampl[digitized == i].mean() for i in range(0, len(mod._xf))]
            ccf1 = self._mod_ccf(mod, ym, y, verbose=False, plot=plot)
            if plot:
                plt.scatter(xmean+xs, ccf1)

            ccf.append(ccf1)

        #plt.plot(mod._xf, y, linewidth=4)
        if weight:
            color = 'r'
        else:
            color = 'g'
        #plt.scatter(xmean+x_shift, ccf/np.max(ccf), c=color)
        try:
            p0 = [np.max(ccf), xmean, 5e-4]
            coeff, var_matrix = curve_fit(gauss, xmean+x_shift, ccf, p0=p0)
            fit = gauss(xmean+x_shift, *coeff)
            ccf_max = coeff[0]
            deltax = coeff[1]-xmean
            deltav = deltax/xmean*aconst.c.to(au.km/au.s).value
            #plt.plot(xmean+x_shift, fit/np.max(fit), c='b')
        except:
            amax = np.argmax(ccf)
            ccf_max = ccf[amax]
            deltax = x_shift[amax]
            deltav = v_shift[amax]
            #plt.scatter(xmean+x_shift[amax], 1)

        if verbose:
            logging.info(("I maximized the data model CCF with a shift of "
                          "%."+str(sd)+"e nm (%."+str(sd)+"e km/s)") \
                          % (deltax, deltav))
        return ccf_max, deltax, deltav


    def _mods_ccf_max(self, vstart, vend, dv, weight):
        systs = self.sess.systs
        for i, m in enum_tqdm(systs._mods_t, len(systs._mods_t),
                              "cookbook_absorbers: Computing CCF"):
            ccf, deltax, deltav = self._mod_ccf_max(m['mod'], vstart, vend, dv,
                                                    weight, verbose=False)
            for i in m['id']:
                w = np.where(systs._t['id']==i)
                systs._t['ccf_deltav'][w] = deltav

        return 0


    def _mods_recreate(self, verbose=True):
        """ Create new system models from a system list """
        spec = self.sess.spec
        spec.t['fit_mask'] = False
        systs = self.sess.systs
        #if len(systs._t)==0: return 0
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        #for i,s in enumerate(systs._t):
        if systs._compressed:
            systs_t = systs._t_uncompressed
        else:
            systs_t = systs._t
        for i,s in enum_tqdm(systs_t, len(systs_t),
                             "cookbook_absorbers: Recreating"):
            systs._id = s['id']
            expr = {}
            for k, v in systs._expr.items():
                #print(k,v)
                if v[0]==systs._id:
                    expr[k] = v[2]
            mod = SystModel(spec, systs, z0=s['z0'], expr=expr)
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
        #print(self.sess.systs._mods_t['mod'])
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
                          None, None, None, None, systs._id])
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
                try:
                    systs._t[iw]['snr'] = np.median(mod._yf*mod._wf)
                except:
                    systs._t[iw]['snr'] = np.nan
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


    def mods_ccf_max(self, vstart=-5, vend=5, dv=0.01, weight=False):
        """ @brief Maximize data/model CCF
        @details Slide the system models around their mean wavelength to
        determine the best data/model CCF and the corresponding shift.
        The wavelength range used for sliding is defined in velocity units.
        @param vstart Range start (km/s with respect to mean wavelength)
        @param vend Range end (km/s with respect to mean wavelength)
        @param dv Range step (km/s)
        @param weight Weight the model by the absolute value of its derivative
        @return 0
        """

        try:
            vstart = float(vstart)
            vend = float(vend)
            dv = float(dv)
            weight = str(weight) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        systs = self.sess.systs
        if 'ccf_deltav' not in systs._t.colnames:
            logging.info("I'm adding column 'ccf_deltav'.")
            systs._t['ccf_deltav'] = np.empty(len(systs._t), dtype=float)

        self._mods_ccf_max(vstart, vend, dv, weight)

        return 0

    def mods_recreate(self):
        """ @brief Recreate the models
        @details Recreate the models from the current system list.
        @return 0
        """

        self._mods_recreate()
        self._spec_update()

        return 0


    def systs_collapse(self):
        """ @brief Collapse the system list
        @details Collapse the list by grouping systems that are modeled together
        @return 0
        """
        self.sess.systs._collapse()
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


    def systs_select(self, z_min=0.0, z_max=10.0, logN_min=10.0, logN_max=18.0,
                     b_min=1.0, b_max=100.0, col=None, col_min=None,
                     col_max=None):
        """ @brief Select systems
        @details Select systems based on their Voigt and fit parameters. A
        logical `and` is applied to all conditions.
        @param z_min Minimum redshift
        @param z_max Maximum redshift
        @param logN_min Minimum (logarithmic) column density
        @param logN_max Maximum (logarithmic) column density
        @param b_min Minimum Doppler broadening
        @param b_max Maximum Doppler broadening
        @param col Other column
        @param col_min Minimum of other column
        @param col_max Maximum of other column
        @return 0
        """

        try:
            z_min = float(z_min)
            z_max = float(z_max)
            logN_min = float(logN_min)
            logN_max = float(logN_max)
            b_min = float(b_min)
            b_max = float(b_max)
            col = None if col in [None, 'None'] else str(col)
            col_min = None if col_min in [None, 'None'] else float(col_min)
            col_max = None if col_max in [None, 'None'] else float(col_max)
        except ValueError:
            logging.error(msg_param_fail)
            return 0


        systs = self.sess.systs

        recompress = False
        if systs._compressed:
            recompress = True
            systs._compress()

        z_sel = np.logical_and(systs._t['z']>z_min, systs._t['z']<z_max)
        logN_sel = np.logical_and(systs._t['logN']>logN_min, systs._t['logN']<logN_max)
        b_sel = np.logical_and(systs._t['b']>b_min, systs._t['b']<b_max)
        cond = np.logical_and(z_sel, np.logical_and(logN_sel, b_sel))

        if col is not None:
            cond = np.logical_and(cond, np.logical_and(systs._t[col]>col_min,
                                                       systs._t[col]<col_max))

        systs._t = systs._t[cond]
        self._mods_recreate()
        self._spec_update()

        if recompress:
            systs._compress()

        return 0


    def systs_sigmav(self):
        """ @brief Estimate position uncertainty
        @details Estimate the uncertainty in the position of systems in velocity
        units.
        @return 0
        """

        spec = self.sess.spec
        systs = self.sess.systs
        lines = self.sess.lines

        if 'fwhm' not in lines._t.colnames:
            logging.error("FWHM of lines is required to compute position "
                          "uncertainty . Please try Recipes > Update lines "
                          "before.")
            return 0

        xpix = np.median(spec._t['xmax']-spec._t['xmin'])

        if 'sigmav' not in systs._t.colnames:
            logging.info("I'm adding column 'sigmav'.")
            systs._t['sigmav'] = at.Column(np.array(np.nan, ndmin=1),
                                           dtype=float)

        for m in systs._mods_t:
            #sel = np.array([np.where(systs._t['id']==id)[0][0] for id in m['id']])
            sel = np.array([], dtype=int)
            for id in m['id']:
                try:
                    sel = np.append(sel, np.where(lines._t['syst_id']==id)[0][0])
                except:
                    pass
            #print(sel)
            #sel = np.array([np.where(lines._t['syst_id']==id)[0][0] for id in m['id']])
            amax = np.argmax(lines._t[sel]['fwhm'])
            fwhm = lines._t[sel]['fwhm'][amax]
            x = lines._t[sel]['x'][amax]
            for i in lines._t[sel]['syst_id']:
                s = np.where(systs._t['id']==i)[0][0]
                systs._t[s]['sigmav'] = (2*np.pi*np.log(2))**(-0.25)/systs._t[s]['snr']\
                                        *np.sqrt(xpix*fwhm)*aconst.c.to(au.km/au.s).value/x

        return 0


    def systs_snr(self):
        """ @brief Estimate SNR of systems
        @details Estimate the signal-to-noise ratio of systems as the median
        flux/flux error ratio in the group interval.
        @return 0
        """

        spec = self.sess.spec
        systs = self.sess.systs

        for m in systs._mods_t:
            for i in m['id']:
                sel = np.where(systs._t['id']==i)[0][0]
                systs._t[sel]['snr'] = np.median(m['mod']._yf*m['mod']._wf)
        return 0


### Advanced

    def cand_find(self, z_start=0, z_end=6, dz=1e-4, resol=resol_def):
        """ @brief Find candidate systems
        @details Cross-match line wavelengths with known transitions to find
        candidate systems.
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param resol Resolution
        @return 0
        """

        try:
            #series = series.replace(';',',')
            z_start = float(z_start)
            z_end = float(z_end)
            dz = float(dz)
            resol = None if resol in [None, 'None'] else float(resol)
        except:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0

        refit_n_temp = dc(self._refit_n)
        max_nfev_temp = dc(self._max_nfev)
        #print(refit_n_temp, max_nfev_temp)

        self._refit_n = 0
        self._max_nfev = 1
        #print(refit_n_temp, max_nfev_temp)

        #for t in np.array(np.meshgrid(trans_d, trans_d)).T.reshape(-1,2):#zip(series_d, series_d):
        trans_arr = np.array(np.meshgrid(trans_d, trans_d)).T.reshape(-1,2)
        for i, t in enum_tqdm(trans_arr, len(trans_arr),
                              "cookbook_absorbers: Finding candidates"):
            p0 = t[0].split('_')[0]
            p1 = t[1].split('_')[0]
            x0 = xem_d[t[0]]
            x1 = xem_d[t[1]]
            z0 = np.min(self.sess.spec.x)/x1-1
            z1 = np.max(self.sess.spec.x)/x0-1
            #print(x0, x1, z0, z1)
            if p0==p1 and x0<x1 and z0>z_start and z1<z_end:
                s = "%s,%s" % (t[0],t[1])
                z_list, logN_list = self.sess.lines._cand_find2(s, z_start, z_end,
                                                        dz)

                s_list = [s]*len(z_list)
                resol_list = [resol]*len(z_list)
                self._systs_add(s_list, z_list, logN_list, resol_list=resol_list, verbose=False)
                #self._systs_cycle()
                self._spec_update()

        self._refit_n = refit_n_temp
        self._max_nfev = max_nfev_temp
        #print(self._refit_n, self._max_nfev)

        return 0


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
