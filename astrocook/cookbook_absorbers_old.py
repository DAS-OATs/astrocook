from .functions import *
from .message import *
from .feat_list import FeatList
from .syst_list import SystList
from .syst_model import SystModel
from .vars import *
from astropy import constants as aconst
from astropy import table as at
from astropy import units as au
from copy import deepcopy as dc
import datetime as dt
import logging
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import argrelmin, argrelmax, find_peaks
from scipy.special import erf, erfc, erfinv
import sys
import time

# from line_profiler import LineProfiler
# from filprofiler.api import profile

prefix = "[INFO] cookbook_absorbers:"

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


class CookbookAbsorbersOld(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        super(CookbookAbsorbersOld, self).__init__()
        self._refit_n = 3
        self._chi2rav_thres = 1e-2
        self._chi2r_thres = np.inf
        self._dlogN_thres = np.inf
        self._max_nfev = max_nfev_def
        self._sel_fit = False

    def _lines_cands_find(self, series, z_start, z_end, dz):
        return self.sess.lines._cands_find(series, z_start, z_end, dz)


    def _logN_guess(self, series, z, b, resol):
        spec = dc(self.sess.spec)
        systs = dc(self.sess.systs)
        ynorm_list = []
        logN_list = np.arange(12, 14, 0.1)
        for logN in logN_list:
            mod = SystModel(systs, z0=z, psf_func = spec.psf_gauss)
            mod._new_voigt(spec, series, z, logN, b, resol,
                           defs=self.sess.defs)
            ynorm_list.append(np.min(mod.eval(x=mod._xs, params=mod._pars)))
        self._guess_f = interp1d(ynorm_list, logN_list-0.5, kind='cubic')


    def _feat_ccf(self, xc, y, ym, verbose=True, plot=False):
        ccf = np.dot(ym, y)
        if plot:
            plt.plot(xc, ym)
            plt.scatter(xc, ym)
        if verbose:
            logging.info("The data-model CCF is %2.3f." % ccf)
        return ccf


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


    def _feat_ccf_max(self, xc, yc, dyc, modelc, vstart=-5, vend=5, dv=1e-2,
                     weight=True, verbose=True):
        sd = -1*int(np.floor(np.log10(dv)))-1
        if sd<0: sd=0

        xmin = xc[0]
        xmax = xc[-1]
        xmean = 0.5*(xmin+xmax)
        v_shift = np.arange(vstart, vend, dv)
        x_shift = xmean * v_shift/aconst.c.to(au.km/au.s).value
        xstart = xmean * vstart/aconst.c.to(au.km/au.s).value
        xend = xmean * vend/aconst.c.to(au.km/au.s).value
        dx = xmean * dv/aconst.c.to(au.km/au.s).value

        #print(vstart, vend, xstart, xend)
        #x_osampl = np.arange(xmin+xstart, xmax+xend, dx)
        x_osampl = np.arange(xmin+xstart, xmax+xend, dx)
        #print()
        #print(x_osampl)
        eval_osampl = 1-np.interp(x_osampl, xc, modelc)
        ccf = []

        y = (1-yc)
        w = 1/dyc
        if weight:
            #w = np.abs(np.gradient(eval_osampl))
            #eval_osampl = eval_osampl * w/np.sum(w)*len(w)
            y = y*w

        #print()
        rc = range(len(xc))
        rc = xc
        #y = (1-mod._yf)#*grad/np.sum(grad)
        #xdiff = np.ediff1d(xc, to_end=np.ediff1d(xc[-2:]))/2
        xdiff = np.ediff1d(xc, to_end=np.ediff1d(xc[-1:]), to_begin=np.ediff1d(xc[:2]))/2
        #print(len(xc), len(x_osampl), len(yc), len(y), len(eval_osampl))
        for i, xs in enumerate(x_shift):
            plot = False
            x = x_osampl+xs
            #digitized = np.digitize(x, xc-xdiff)-1
            #ym = [eval_osampl[digitized == j].mean() for j in range(len(xc))]
            ym = np.interp(xc, x, eval_osampl)
            ccf1 = self._feat_ccf(xc, y, ym, verbose=False, plot=plot)
            if plot:
                plt.scatter(xmean+xs, ccf1/500)

            ccf.append(ccf1)

        #plt.plot(rc, y, linewidth=2, c='b')
        #plt.scatter(rc, 1-modelc, linewidth=2, c='lightblue')
        #plt.plot(x_osampl, eval_osampl, linewidth=1, c='green')
        if weight:
            color = 'r'
        else:
            color = 'black'
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
            #print(ccf)
            amax = np.argmax(ccf)
            ccf_max = ccf[amax]
            deltax = x_shift[amax]
            deltav = v_shift[amax]
            #plt.scatter(xmean+x_shift[amax], 1)
        #print(len(xc), len(xdiff))
        #print(np.any(np.ediff1d(xc-xdiff)<=0))
        #print(np.ediff1d(xc-xdiff))
        xshift = x_osampl+deltax
        digitized = np.digitize(xshift, xc-xdiff)-1
        yshift = 1-np.array([eval_osampl[digitized == j].mean() for j in range(len(xc))])

        if deltav < vstart or deltav > vend:
            ccf_max = np.nan
            deltax = np.nan
            deltav = np.nan
        #print()
        #print(deltax,deltav)
        if verbose:
            logging.info(("I maximized the data model CCF with a shift of "
                          "%."+str(sd)+"e nm (%."+str(sd)+"e km/s)") \
                          % (deltax, deltav))
        if plot: plt.show()
        return ccf_max, deltax, deltav, yshift


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


    def _feats_ccf_max(self, vstart, vend, dv, weight, xcol='x', ycol='y',
                       dycol='dy', contcol='cont', modelcol='model', thr=1e-1,
                       update_modelcol=False, height=1e-1, prominence=1e-1):
        weight = str(weight) == 'True'
        spec = self.sess.spec
        systs = self.sess.systs
        feats = np.hstack(([0], systs._bounds, [-1]))

        deltav_arr = np.array([])
        xmean_arr = np.array([])
        if update_modelcol:
            spec._t[modelcol+'_shift'] = spec._t[modelcol]
        xcf = []
        ycf = []
        contcf = []
        feats = FeatList() #xunit=spec._xunit, yunit=spec._yunit)
        feats._maxs_from_spec(spec, height, prominence)
        #print(feats._maxs)
        for i, f in enum_tqdm(feats._maxs[:-1], len(feats._maxs)-1,
                              "cookbook_absorbers: Computing CCF for features"):
            fe = feats._maxs[i+1]
            sel = np.s_[f:fe]
            cut = np.where(spec._t[modelcol][sel]/spec._t[contcol][sel]<1-thr)
            contc = spec._t[contcol][sel][cut]
            xc = spec._t[xcol][sel][cut]
            yc = spec._t[ycol][sel][cut]/contc
            dyc = spec._t[dycol][sel][cut]/contc
            modelc = spec._t[modelcol][sel][cut]/contc
            if len(cut[0])>0:
                feats._add(spec._t[sel][cut], systs)
                ccf, deltax, deltav, modelshift = self._feat_ccf_max(xc, yc, dyc, modelc,
                                                         vstart, vend, dv,
                                                         weight, verbose=False)
                feats._l[-1]._ccf_compute(np.mean(xc), deltav)
            #else:
                #print('xc len 0')
            #    deltav = 0.0
            #    modelshift = modelc
                xmean_arr = np.append(xmean_arr, np.mean(xc))
                deltav_arr = np.append(deltav_arr, deltav)
                """
            for i in m['id']:
                w = np.where(systs._t['id']==i)
                systs._t['ccf_deltav'][w] = deltav
                """
                if update_modelcol:
                    spec._t[modelcol+'_shift'][sel][cut] = modelshift*spec._t[contcol][sel][cut]

        empty = [np.nan]*len(xcf)
        feats._table_update()
        self.sess.feats = feats
        #print(len(deltav_arr), len(feats._l))
        with open(self.sess.name+'_deltav.npy', 'wb') as f:
            np.save(f, xmean_arr)
            np.save(f, deltav_arr)
        return 0

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

    def _mods_recreate(self, **kwargs):
        return self._mods_recreate2(**kwargs)

    def _mods_recreate1(self, verbose=True):
        """ Create new system models from a system list """
        spec = self.sess.spec
        spec.t['fit_mask'] = False
        systs = self.sess.systs
        #if len(systs._t)==0: return 0
        systs._mods_t.remove_rows(range(len(systs._mods_t)))
        #for i,s in enumerate(systs._t):
        compressed = False
        if systs is not None and systs._compressed:
            systs_t = systs._t_uncompressed
        else:
            systs_t = systs._t
        for i,s in enum_tqdm(systs_t, len(systs_t),
                             "cookbook_absorbers: Recreating"):
            systs._id = s['id']
            vars = {}
            constr = {}
            #print(systs._constr)
            for k, v in systs._constr.items():
                if v[0]==systs._id:
                    if v[2]!=None:
                        constr[k] = v[2]
                    else:
                        vars[k.split('_')[-1]+'_vary'] = False
            mod = SystModel(systs, z0=s['z0'], vars=vars, constr=constr, psf_func = spec.psf_gauss)
            mod._new_voigt(spec, series=s['series'], z=s['z'], logN=s['logN'],
                           b=s['b'], resol=s['resol'],
                           defs=self.sess.defs)
            self._mods_update(mod)
        mods_n = len(self.sess.systs._mods_t)
        if verbose:
            logging.info("I've recreated %i model%s." \
                         % (mods_n, '' if mods_n==1 else 's'))
        return 0

    def _mods_recreate2(self, only_constr=False, mod_new=None, rem_id=None,
                        fast=False, verbose=True):
        """ Create new system models from a system list """
        spec = self.sess.spec
        spec.t['fit_mask'] = False
        systs = self.sess.systs

        # When constraints have been added
        if only_constr:
            mod_sel = np.array([], dtype=int)
            mod_w = np.array([], dtype=int)
            for k,v in systs._constr.items():
                if v[2]!=None:
                    if v[2]=='':
                        mod_sel = np.append(mod_sel, v[0])
                    else:
                        mod_sel = np.append(mod_sel,
                                            [int(v[2].split('_')[-2]),v[0]])
            mod_sel = np.ravel(mod_sel)

            for i in range(2):
                for id in mod_sel:
                    mod_w = np.append(mod_w, np.where([id in mod_id for mod_id in systs._mods_t['id']])[0])
                mod_w = np.unique(mod_w)
                mod_sel = np.array([], dtype=int)
                for w in mod_w:
                    mod_sel = np.append(mod_sel, np.array([systs._mods_t['id'][w]]))

        # When a model has been added
        elif mod_new is not None:
            mod_sel = np.array([], dtype=int)
            mod_w = np.array([], dtype=int)
            ys = mod_new._ys
            for s in systs._mods_t:
                mod = s['mod']
                ys_s = mod._ys
                ymax = np.maximum(ys, ys_s)
                thres = 1e-2
                y_cond = np.amin(ymax)<1-thres or np.amin(ymax)==np.amin(ys)
                pars_cond = False
                for p,v in mod_new._constr.items():
                    for mod_p,mod_v in mod._pars.items():
                        pars_cond = pars_cond or v==mod_p
                if y_cond or pars_cond:
                    for id in s['id']:
                        mod_w = np.append(mod_w,
                            np.where([id in mod_id for mod_id
                                      in systs._mods_t['id']])[0])
                    mod_w = np.unique(mod_w)
                    mod_sel = np.array([], dtype=int)
                    for w in mod_w:
                        mod_sel = np.append(mod_sel, np.array([systs._mods_t['id'][w]]))

        #When a system has been removed
        elif rem_id is not None:
            mod_sel = np.array([], dtype=int)
            t_id = systs._t['id']
            mods_t_id = systs._mods_t['id']
            for r in rem_id:
                mod_w = [r in m for m in mods_t_id]
            mod_sel = np.append(mod_sel, mods_t_id[mod_w][0])


        else:
            mod_w = range(len(systs._mods_t))
            mod_sel = np.array(systs._t['id'])

        compressed = False
        if systs is not None and systs._compressed:
            systs_t = systs._t_uncompressed
        else:
            systs_t = systs._t

        # Collect existing specifications for logN_tot
        N_tot_specs_d = {}
        if hasattr(systs, '_N_tot_specs'):
            N_tot_specs_d = systs._N_tot_specs

        if not fast:
            systs._mods_t.remove_rows(mod_w)

            systs_t.sort('id')
            wrong_id = []
            corr_id = []
            tt = time.time()
            for i,s in enum_tqdm(systs_t, len(systs_t), #len(mod_sel),
                                 "cookbook_absorbers: Recreating"):
            #for i,m in enum_tqdm(mod_sel, len(mod_sel),
            #                     "cookbook_absorbers: Recreating"):
                systs._id = s['id']
                #systs._id = m
                if systs._id in mod_sel:
                #for s in systs_t[np.where(systs_t['id']==m)[0]]:
                    vars = {}
                    constr = {}

                    for k, v in systs._constr.items():
                        if v[0]==systs._id:
                            try:
                                id_check = v[2].split('_')[-2]
                            except:
                                id_check = ''
                            ids = [str(id) for id in systs_t['id']]+['']

                            if v[2]!=None and id_check in ids:
                                constr[k] = v[2]
                            else:
                                vars[k.split('_')[-1]+'_vary'] = False
                    mod = SystModel(systs, z0=s['z0'], vars=vars, constr=constr, psf_func = spec.psf_gauss)
                    if any([mod._id in i for i in systs._mods_t['id']]):
                        wrong_id.append(mod._id)
                        corr_id.append(np.max(systs_t['id'])+1)
                        mod._id = np.max(systs_t['id'])+1

                    # Implement existing specifications for logN_tot
                    if mod._id in N_tot_specs_d:
                        N_tot = True
                        N_tot_specs = N_tot_specs_d[mod._id]
                    else:
                        N_tot = False
                        N_tot_specs = (None, None, None)
                    mod._new_voigt(spec, series=s['series'], z=s['z'], logN=s['logN'],
                                   b=s['b'], resol=s['resol'],
                                   defs=self.sess.defs,
                                   N_tot=N_tot, N_tot_specs=N_tot_specs)
                    self._mods_update(mod)

                else:
                    systs._id = np.max(systs._t['id'])+1

            #systs._id = np.max(systs._t['id'])+1

            for w, c in zip(wrong_id, corr_id):
                logging.warning("System %i had a duplicated id! I changed it "
                                "to %i." % (w, c))
        else:
            systs._id = np.max(systs._t['id'])+1

        spec.t['fit_mask'] = False
        active_c = 0
        for i, m in enumerate(systs._mods_t):
            mod = m['mod']
            try:
                active = mod._active
            except:
                active = True
            if active:
                active_c += 1
                c = []
                t = 1e-2
                j = 0
                while len(c)==0 and j<1000:
                    try:
                        c = np.where(mod._yl<1-t)[0]
                    except:
                        c = np.where(mod._ys<1-t)[0]
                    t = t*0.5
                    #print(i, j, len(c))
                    j += 1

                spec.t['fit_mask'][c] = True


        systs_t.sort(['z','id'])
        systs_n = len(systs._t)
        mods_n = len(systs._mods_t)
        if verbose:
            logging.info("I've recreated %i model%s (including %i system%s)." \
                         % (mods_n, '' if mods_n==1 else 's',
                            systs_n, '' if systs_n==1 else 's'))
            if active_c < mods_n:
                logging.info("Only %i model%s %s active and eligible for "
                             "fitting." % (active_c, '' if active_c==1 else 's',
                             'is' if active_c==1 else 'are'))


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
        systs._mods_t.sort('id')
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
        #print(cont[s],y[s],model[s])
        systs._bounds = argrelmax(model[s]/cont[s])[0]
        try:
            deabs[s] = cont[s] + y[s] - model[s]
        except:
            deabs[s] = np.array(cont[s]) + np.array(y[s]) - np.array(model[s])
        return 0


    def _syst_add(self, series, z, logN, b, resol, verbose=True):
        systs = self.sess.systs
        spec = self.sess.spec
        #print(systs._t['series'][systs._t['z0']==z])
        if z in systs._t['z0'] \
            and series in systs._t['series'][systs._t['z0']==z]:
            if verbose:
                logging.warning("Redshift %2.4f already exists. Choose another "
                                "one." % z)
            return None
        systs._t.add_row(['voigt', series, z, z, None, logN, None, b,
                          None, 0.0, None, resol, None, None, systs._id])
        #systs._id = np.max(systs._t['id'])+1
        from .syst_model import SystModel
        mod = SystModel(systs, z0=z, psf_func = spec.psf_gauss)
        #print(self.sess.defs.dict['voigt'])
        mod._new_voigt(spec, series, z, logN, b, resol, defs=self.sess.defs)

        # When a single system is added, it is stored only on the model table
        self._mods_update(mod, incr=False)
        return mod


    def _syst_fit(self, mod, verbose=True):
        if self._max_nfev > 0:
            try:
                if 'max_nfev' in self.sess.defs.dict['fit']:
                    fit_kws = self.sess.defs.dict['fit']
                else:
                    fit_kws = self.sess.defs.dict['fit']
                    fit_kws['max_nfev'] = self._max_nfev
            except:
                fit_kws = {'max_nfev': self._max_nfev}
            #print(fit_kws)
            #print(self.sess.defs.dict['fit'])
            if mod._active:
                frozen = mod._fit(fit_kws=fit_kws)
            else:
                frozen = 1
            #mod._pars.pretty_print()
            if verbose and frozen:
                if not mod._active:
                    logging.info("I've not fitted 1 model at redshift %2.4f "
                                 "because it is not active." % mod._z0)
                else:
                    logging.info("I've not fitted 1 model at redshift %2.4f "
                                 "because all the parameters are frozen."
                                 % mod._z0)
            elif verbose:
                logging.info("I've fitted 1 model at redshift %2.4f." \
                             % mod._z0)
        else:
            frozen = 1
            logging.info("I'm not fitting the system because you choose "
                         "max_nfev=0.")

        # When a single system is fitted, it is stored also the system table
        if not frozen: self._systs_update(mod)
        return frozen


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
                   resol_list=None, k_list=None, verbose=True):
        if logN_list is None: logN_list = [None]*len(series_list)
        if b_list is None: b_list = [None]*len(series_list)
        if resol_list is None: resol_list = [None]*len(series_list)
        if k_list is None: k_list = [None]*len(series_list)
        systs_n = 0
        id_list = []
        for i, (series, z, logN, b, resol, k) \
            in enum_tqdm(zip(series_list, z_list, logN_list, b_list, resol_list, k_list),
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
                id_list.append(mod._id)
                if k is not None:
                    self.sess.systs._constr['lines_voigt_%i_z' % mod._id] \
                        = (mod._id, 'z', k)


        # Improve
        mods_t = self.sess.systs._mods_t

        if verbose and systs_n>0:
            if len(np.unique(series_list))==1:
                logging.info("I've added %i %s system%s in %i model%s." \
                             % (systs_n, series_list[0],
                             '' if systs_n==1 else 's',
                             len(mods_t), msg_z_range(z_list)))
            else:
                logging.info("I've added %i system%s in %i model%s." \
                             % (systs_n, '' if systs_n==1 else 's',
                             len(mods_t), msg_z_range(z_list)))
        return id_list


    def _systs_compress(self):
        self.sess.systs._compress()
        return 0


    def _systs_cycle(self, mod=None, verbose=True, recreate=True):
        chi2rav = np.inf
        chi2rav_old = 0
        chi2r_list, z_list = [], []
        for i,_ in enum_tqdm(range(self._refit_n), self._refit_n,
                              'cookbook_absorbers: Cycling'):
            if chi2rav > self._chi2rav_thres and chi2rav != chi2rav_old:
                if chi2rav < np.inf: chi2rav_old = chi2rav
                fit_list, chi2r_list, z_list = self._systs_fit(verbose=False)
                if i > 1 and len(chi2r_list)==len(chi2r_list_old):
                    chi2rav = np.mean(np.abs(np.array(chi2r_list)\
                                             -np.array(chi2r_list_old)))
                chi2r_list_old = chi2r_list
                self._systs_reject(mod=mod, verbose=verbose)
                #self._mods_recreate(verbose=False)
                self._mods_recreate(mod_new=mod, verbose=verbose)
            #print(chi2rav, chi2rav_old)
        fit_list, chi2r_list, z_list = self._systs_fit(verbose=False)

        self._systs_reject(mod=mod, verbose=verbose)
        if recreate: self._mods_recreate(mod_new=mod, verbose=verbose)
        if verbose and z_list != []:
            logging.info("I've fitted %i model%s." \
                         % (np.sum(fit_list), msg_z_range(z_list)))
            if chi2rav < np.inf:
                logging.info("Average chi2r variation after last cycle: %2.4e."\
                             % chi2rav)


    def _systs_fit(self, verbose=True):
        systs = self.sess.systs
        mods_t = systs._mods_t
        z_list = []
        chi2r_list = []
        if self._max_nfev > 0:
            fit_list = []
            for i,m in enumerate(mods_t):
                if self._sel_fit:
                    dz = [systs._t['dz'][np.where(systs._t['id']==id)[0][0]] \
                          for id in m['id']]
                    fit_list.append(np.isnan(dz).any())
                else:
                    fit_list.append(True)

            for i,m in enum_tqdm(mods_t, np.sum(fit_list),
                                 "cookbook_absorbers: Fitting"):
            #for i,m in enumerate(mods_t):
                """
                if self._sel_fit:
                    dz = [systs._t['dz'][np.where(systs._t['id']==id)[0][0]] \
                          for id in m['id']]
                    fit = np.isnan(dz).any()
                else:
                    fit = True
                """
                if fit_list[i]:
                    z_list.append(m['z0'])
                    frozen = self._syst_fit(m['mod'], verbose=False)
                    if frozen:
                        fit_list[i] = False
                    else:
                        chi2r_list.append(m['mod']._chi2r)
            if verbose:
                logging.info("I've fitted %i model%s." \
                             % (np.sum(fit_list), msg_z_range(z_list)))
        else:
            fit_list = []
            for i,m in enumerate(mods_t):
                z_list.append(m['z0'])
                fit_list.append(False)
            if verbose:
                logging.info("I've not fitted any model because you choose "
                             "max_nfev=0.")
        return fit_list, chi2r_list, z_list


    def _systs_guess(self, series_list, z_list):
        logN_list = np.array([])
        for series, z in zip(series_list, z_list):
            logN_list = np.append(logN_list, self._syst_guess(series, z))
        return logN_list


    def _systs_unify(self, dz=1e-5):
        systs = dc(self.sess.systs)
        merged = np.array([], dtype=int)
        for syst in systs._t:
            z = syst['z']
            pref = syst['series'].split('_')[0]
            systs_pref = np.array([s.split('_')[0] for s in systs._t['series']])
            w = np.where(np.logical_and(systs._t['z']>z-dz,
                                        systs._t['z']<z+dz))[0]
            sel = systs_pref[w]==pref
            if np.sum(sel)>1:
                #print(w[sel])
                #print(systs._t['series'][w[sel][0]])
                self.sess.systs._t['series'][w[sel][0]] = \
                    ','.join(np.unique(systs._t['series'][w[sel]]))
                merged = np.append(merged, w[sel][1:])

        if len(merged)>0:
            #print(merged)
            self.sess.systs._t.remove_rows(merged)
            self._mods_recreate()

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

    def _systs_reject(self, mod=None, verbose=True):
        systs = self.sess.systs
        chi2r_cond = systs._t['chi2r'] > self._chi2r_thres
        """
        relerr_cond = np.logical_or(np.logical_or(
            systs._t['dz'] > dlogN_thres*systs._t['z'],
            systs._t['dlogN'] > dlogN_thres*systs._t['logN']),
            systs._t['db'] > dlogN_thres*systs._t['b'])
        """
        relerr_cond = systs._t['dlogN'] > self._dlogN_thres
        cond = np.logical_or(chi2r_cond, relerr_cond)
        if mod is not None:
            id_cond = np.array([mod._id==ids for ids in systs._t['id']])
            chi2r_cond = np.logical_and(chi2r_cond, id_cond)
            relerr_cond = np.logical_and(relerr_cond, id_cond)
            cond = np.logical_and(cond, id_cond)
        rem = np.where(cond)[0]
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
                logging.info("I've rejected %i badly fitting system%s (%i with a "\
                             "reduced chi2 above %2.2f, %i with relative errors "\
                             "above %2.2f)."
                             % (len(rem), '' if len(rem)==1 else 's',
                                np.sum(chi2r_cond), self._chi2r_thres,
                                np.sum(relerr_cond), self._dlogN_thres))
        return 0 #refit_id


    def _systs_remove(self, rem):#, refit_id):
        systs = self.sess.systs
        t_id = systs._t['id']
        mods_t_id = systs._mods_t['id']
        rem_id = [t_id[r] for r in rem]
        for i, r in enum_tqdm(rem, len(rem), "cookbook_absorbers: Removing"):
            sel = [t_id[r] in m for m in mods_t_id]
            k_del = []
            for k, v in systs._constr.items():
                try:
                    id_check = v[2].split('_')[-2]
                except:
                    id_check = ''
                if str(t_id[r]) in [v[0], id_check]:
                    k_del.append(k)

        systs._t.remove_rows(rem)
        for k in k_del:
            del systs._constr[k]

        self._mods_recreate(rem_id=rem_id)
        return 0


    def _systs_update(self, mod, incr=True):
        systs = self.sess.systs
        #print(systs._mods_t['mod'])
        modw = np.where(mod == systs._mods_t['mod'])[0][0]
        ids = systs._mods_t['id'][modw]
        for i in ids:
            try:
                iw = np.where(systs._t['id']==i)[0][0]
                pref = 'lines_voigt_'+str(i)
                systs._t[iw]['z'] = mod._pars[pref+'_z'].value
                systs._t[iw]['dz'] = mod._pars[pref+'_z'].stderr
                if pref+'_N_tot' in mod._pars:
                    logN = np.log10(mod._pars[pref+'_N_tot'].value\
                                    -mod._pars[pref+'_N_other'].value)
                    if np.isnan(logN):
                        logN = 10
                    systs._t[iw]['logN'] = logN
                    systs._t[iw]['dlogN'] = np.nan
                else:
                    systs._t[iw]['logN'] = mod._pars[pref+'_logN'].value
                    systs._t[iw]['dlogN'] = mod._pars[pref+'_logN'].stderr
                systs._t[iw]['b'] = mod._pars[pref+'_b'].value
                systs._t[iw]['db'] = mod._pars[pref+'_b'].stderr
                systs._t[iw]['btur'] = mod._pars[pref+'_btur'].value
                systs._t[iw]['dbtur'] = mod._pars[pref+'_btur'].stderr
                try:
                    systs._t[iw]['resol'] = mod._pars['psf_gauss_%i_resol' % i].value
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


    def series_replace(self, series_from='Ly_a', series_to='Ly_a,Ly_b'):
        """ @brief Replace series
        @details Replace all given series in a table with a different series
        @param series_from Replace series
        @param series_to with
        @return 0
        """

        try:
            series_from = str(series_from)
            series_to = str(series_to)
        except:
            logging.error(msg_param_fail)
            return 0

        systs = self.sess.systs
        w = systs._t['series']==series_from
        systs._t['series'][w] = series_to
        self._mods_recreate()
        self._spec_update()


    def comp_extract(self, num=1):
        """ @brief Extract systems based on components
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


    def feats(self, thres=1e-2, height=1e-1, prominence=1e-2):
        """ @brief Organize systems into absorption features
        @details Absorption features are portions of the model including one or
        more systems, and limited by local maxima in the model.
        @param thres Threshold for cutting the model when it gets close to continuum, normalized to continuum
        @param height Minimum height for the local maxima to count as boundaries, normalized to continuum
        @param prominence Minimum prominence for the local maxima to count as boundaries, normalized to continuum
        @return 0
        """

        try:
            thres = float(thres)
            height = float(height)
            prominence = float(prominence)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.feats = FeatList()
        self.sess.feats.create(self.sess.spec, self.sess.systs, thres, height,
                               prominence)
        logging.info("I've extracted %i features from the system models."\
                     % len(self.sess.feats._l))

        return 0


    def feats_z_lock(self):
        """ @brief Lock redshifts by absorption features
        @details Lock the difference of systems redshift by absorption features.
        @return 0
        """

        self.sess.feats._z_lock()
        logging.info("I've locked redshifts by absorption features.")

        return 0


    def feats_logN_tot(self, set_specs=True, sel=None):
        """ @brief Prepare features to fit total column density
        @details The last system in each feature is prepared to be modeled as
        a function of the total column density of the feature itself, to obtain
        realistic errors on the total column density.
        @param set_pars Set the parameters to model the systems (if False, total column densities will be only extracted from the models)
        @param sel Feature selected to be modeled in this way (None to select all)
        @return 0
        """

        self.sess.feats._logN_tot(self.sess.systs, set_specs, sel)
        if set_specs:
            self._mods_recreate()
            for m in self.sess.systs._mods_t:
                #m._pars.pretty_print()
                self._systs_update(m['mod'])
                for i in m['id']:
                    self.sess.systs._d[i]._mod = m['mod']

            #for s in self.sess.systs._d:
            #    self.sess.systs._d[s]._mod._pars.pretty_print()
            self.sess.feats._systs_update(self.sess.systs)
        #for f in self.sess.feats._l:
        #    s = max(f._systs.keys())
        logging.info("I've prepared features to fit total column density.")

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

        return

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


    def chunk_fit(self, #xmin, xmax,
                  chunks, refit_n=0, chi2rav_thres=1e-2,
                  max_nfev=max_nfev_def, recreate=True):
        """ @brief Fit systems in a spectrum chunk
        @details Fit one or more system in a spectrum chunk, freezing the
        components of all other systems.
        @param chunks Wavelength range of the chunks (nm), e.g. 500-501;502-503
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @return 0
        """
        try:
            #xmin = float(xmin) * au.nm
            #xmax = float(xmax) * au.nm
            chunks = chunk_parse(chunks)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        start = dt.datetime.now()

        #chunks = chunk_parse(chunks)

        t = self.sess.systs._t
        x = np.array([np.array([to_x(z, trans).value for trans in trans_parse(s)])
                      for (z,s) in t['z', 'series']])
        w_all = np.array([], dtype=int)
        for c in chunks:
            w = np.where([np.any(np.logical_and(xi>c[0],xi<c[1])) for xi in x])[0]
            w_all = np.append(w_all, w)
        ids = t['id'][w_all]
        self.sess.systs._freeze_pars(exclude=ids)

        # Select model
        mods_t = self.sess.systs._mods_t

        self.systs_fit(refit_n, recreate=recreate)
        end = dt.datetime.now()

        chunks_s = '['
        chunks_s += '] nm, ['.join(['%3.3f-%3.3f' % (c[0], c[1]) for c in chunks])
        chunks_s += '] nm'
        logging.info("I've fitted %i system%s in %s in %3.3f seconds." \
                     % (len(ids), 's' if len(ids)>1 else '', chunks_s,
                        (end-start).total_seconds()))

        self.sess.systs._unfreeze_pars(exclude=ids)
        self._spec_update()

        return 0



    def syst_fit(self, ids=[1], refit_n=0, chi2rav_thres=1e-2,
                 max_nfev=max_nfev_def):
        """ @brief Fit individual systems
        @details Fit one or more system, freezing the components of all other systems.
        @param id System id
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @return 0
        """
        try:
            ids = [int(i) for i in ids[1:-1].split(',')]
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.systs._freeze_pars(exclude=ids)

        # Select model
        mods_t = self.sess.systs._mods_t

        for i,id in enum_tqdm(ids, len(ids), 'cookbook_absorbers: Fitting'):
            for m in mods_t:
                if id in m['id']:
                    mod = m['mod']
            self.systs_fit(refit_n, _mod=mod)
        logging.info("I've fitted %i system%s." \
                     % (len(ids), 's' if len(ids)>1 else ''))

        self.sess.systs._unfreeze_pars(exclude=ids)
        self._spec_update()

        return 0

    def group_extract(self, id=1):
        """ @brief Extract a group of systems
        @details Extract all system grouped to the selected system and remove
        all the others.
        @param id System id
        @return 0
        """

        try:
            id = int(id)
        except:
            logging.error(msg_param_fail)
            return 0

        out = self.sess.systs
        t_sel = []
        mods_t_sel = []
        for i, m in enumerate(out._mods_t):
            if id not in m['id']:
                for mid in m['id']:
                    t_sel.append(np.where(out._t['id'] == mid)[0][0])
                mods_t_sel.append(i)
        out._t.remove_rows(t_sel)
        out._mods_t.remove_rows(mods_t_sel)

        return 0


    def group_fit(self, id=1, refit_n=0, chi2rav_thres=1e-2,
                  max_nfev=max_nfev_def):
        """ @brief Fit a group of systems
        @details Fit together all systems that are grouped to the selected system.
        @param id System id
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @return 0
        """
        try:
            id = int(id)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        # Select model
        mods_t = self.sess.systs._mods_t
        for m in mods_t:
            if id in m['id']:
                mod = m['mod']

        # Select system
        t = self.sess.systs._t
        for s in t:
            if id == s['id']:
                s['dz'] = np.nan
        self._sel_fit = True

        self._systs_cycle(mod=mod)
        self._spec_update()

        return 0


    def syst_fit_old(self, num=1, refit_n=0, chi2rav_thres=1e-2,
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
            num = int(num)-1
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
        except:
            logging.error(msg_param_fail)
            return 0

        mods_t = self.sess.systs._mods_t
        mod = mods_t['mod'][num in mods_t['id']]
        #self._syst_fit(mod)
        self._systs_cycle(mod)
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


    def systs_fit(self, refit_n=3, chi2rav_thres=1e-2, max_nfev=max_nfev_def,
                  sel_fit=False, _mod=None, recreate=True):
        """ @brief Fit systems
        @details Fit all Voigt model from a list of systems.
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @param sel_fit Selective fit (only new systems will be fitted)
        @return 0
        """

        try:
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = int(max_nfev)
            self._sel_fit = str(sel_fit) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        #self._systs_fit()
        self._systs_cycle(mod=_mod, verbose=False, recreate=recreate)
        self._spec_update()

        return 0

    def systs_supersede(self, dv=5, series='Ly-a'):
        """ @brief Supersede systems
        @details Enforce rules for systems to supersede other systems in case of
        superposition. Two systems are superposed when the difference between
        their positions, expressed as a velocity, is below a given threshold.
        @param dv Velocity threshold (km/s)
        @param series Series to be superseded
        @return 0
        """

        try:
            dv = float(dv)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        systs = self.sess.systs


        q = [trans_parse(t) for t in systs._t['series']]
        #print(systs._t['z'], s)
        r = np.array([], dtype=int)
        #s = np.array([])
        x = np.array([])
        for i, (zi, si) in enumerate(zip(systs._t['z'], systs._t['series'])):
            r = np.append(r, [i for t in trans_parse(si)])
            #s = np.append(s, trans_parse(si))
            x = np.append(x, [to_x(zi, t).value for t in trans_parse(si)])
        argsort = np.argsort(x)
        v = aconst.c.to(au.km/au.s).value*x/x[argsort][0]
        where = np.where(np.ediff1d(v[argsort])<dv)

        rem = []
        for w in where[0]:
            check_1 = systs._t['series'][r[argsort][w]]!= series
            check_2 = systs._t['series'][r[argsort][1:][w]]!= series
            if check_1 and not check_2:
                rem.append(r[argsort][1:][w])
            if check_2 and not check_1:
                rem.append(r[argsort][w])

        self._systs_remove(rem)
        self._mods_recreate()

        return 0

    def _feats_select_isolated(self, thres=1e-1):
        if not hasattr(self.sess, 'feats'):
            self.feats()
        self.sess.feats._select_isolate(thres)

        return 0


    def _feats_select(self, z_min=0.0, z_max=10.0, logN_min=10.0, logN_max=22.0,
                      b_min=1.0, b_max=100.0, col=None, col_min=None,
                      col_max=None):

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


        spec = self.sess.spec
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        z_sel = np.logical_and(systs._t['z']>z_min, systs._t['z']<z_max)
        logN_sel = np.logical_and(systs._t['logN']>logN_min, systs._t['logN']<logN_max)
        b_sel = np.logical_and(systs._t['b']>b_min, systs._t['b']<b_max)
        cond = np.logical_and(z_sel, np.logical_and(logN_sel, b_sel))

        if col is not None:
            cond = np.logical_and(cond, np.logical_and(systs._t[col]>col_min,
                                                       systs._t[col]<col_max))

        xs = np.array([])
        for s in systs._t[cond]:
            for si in trans_parse(s['series']):
                xs = np.append(xs, to_x(s['z'],si).value)

        feats = np.hstack(([0], systs._bounds, [-1]))
        sel = np.histogram(xs, spec._t['x'][feats])[0]>0
        with open(self.sess.name+'_feats_sel.npy', 'wb') as f:
            np.save(f, sel)

        return 0


    def systs_select(self, series='any', z_min=0.0, z_max=10.0, logN_min=10.0,
                     logN_max=22.0, b_min=1.0, b_max=100.0, col=None,
                     col_min=None, col_max=None):
        """ @brief Select systems
        @details Select systems based on their Voigt and fit parameters. A
        logical `and` is applied to all conditions.
        @param series Series
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
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        if series != 'any':
            series_sel = [np.any([t in trans_parse(series)
                                  for t in trans_parse(s['series'])])
                          for s in systs._t]
        else:
            series_sel = np.ones(len(systs._t))
        z_sel = np.logical_and(systs._t['z']>z_min, systs._t['z']<z_max)
        logN_sel = np.logical_and(systs._t['logN']>logN_min, systs._t['logN']<logN_max)
        b_sel = np.logical_and(systs._t['b']>b_min, systs._t['b']<b_max)
        cond = np.logical_and(np.logical_and(series_sel, z_sel), np.logical_and(logN_sel, b_sel))

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

    def cands_find(self, series='all', z_start=0, z_end=6, dz=1e-4,
                  resol=resol_def, avoid_systs=True, append=True):
        """ @brief Find candidate systems
        @details Cross-match line wavelengths with known transitions to find
        candidate systems.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param resol Resolution
        @param avoid_systs Avoid finding candidates over systems already
        detected
        @param append Append systems to existing system list
        @return 0
        """

        try:
            #series = series.replace(';',',')
            #series = None if series in [None, 'None'] else str(series)
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            resol = None if resol in [None, 'None'] else float(resol)
            avoid_systs = str(avoid_systs) == 'True'
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check_old(self.sess.spec, resol)
        if not check: return 0

        refit_n_temp = dc(self._refit_n)
        max_nfev_temp = dc(self._max_nfev)
        #print(refit_n_temp, max_nfev_temp)

        self._refit_n = 0
        self._max_nfev = 1
        #print(refit_n_temp, max_nfev_temp)

        #for t in np.array(np.meshgrid(trans_d, trans_d)).T.reshape(-1,2):#zip(series_d, series_d):
        #if short:
        #    t_d = trans_d_short
        if series != 'all':
            t_d = trans_parse(series)
        else:
            t_d = trans_d
        trans_arr = np.array(np.meshgrid(t_d, t_d)).T.reshape(-1,2)
        z_list = np.array([])
        logN_list = np.array([])
        s_list = np.array([])
        resol_list = np.array([])
        count = 0
        for i, t in enum_tqdm(trans_arr, len(trans_arr),
                              "cookbook_absorbers: Finding candidates"):
            p0 = t[0].split('_')[0]
            p1 = t[1].split('_')[0]
            x0 = xem_d[t[0]]
            x1 = xem_d[t[1]]
            z0 = np.min(self.sess.spec.x)/x1-1
            z1 = np.max(self.sess.spec.x)/x0-1
            z_start = max(z0, z_start)
            z_end = min(z1, z_end)
            if p0==p1 and x0<x1:# and (z0>z_start or z1<z_end):
                s = "%s,%s" % (t[0],t[1])
                #print(s)
                z_l, logN_l, _ = self.sess.lines._cands_find2(s, z_start, z_end, dz)

                self._systs_prepare(append)

                if avoid_systs and 'model' in self.sess.spec.t.colnames:
                    for zi in z_l:
                        for si in trans_parse(s):
                            xi = to_x(zi,si)
                            wi = np.abs(self.sess.spec.x - xi).argmin()
                            mi = self.sess.spec.t['model'][wi]
                            #add = bool(add and mi>1-1e-4)
                            if mi<1-1e-4:
                                wi = np.where(z_l!=zi)[0]
                                z_l = np.array(z_l)[wi]
                                logN_l = np.array(logN_l)[wi]

                add = len(z_l)>0
                s_l = [s]*len(z_l)
                resol_l = [resol]*len(z_l)

                if add:
                    #z_list = np.append(z_list, z_l)
                    #logN_list = np.append(logN_list, logN_l)
                    #s_list = np.append(s_list, s_l)
                    #resol_list = np.append(resol_list, resol_l)
                    count += len(z_l)
                    self._systs_add(s_l, z_l, logN_l, resol_list=resol_l, verbose=False)
                self._spec_update()
        #self._systs_prepare(append)
        #self._systs_add(s_list, z_list, logN_list, resol_list=resol_list, verbose=False)
        #self._spec_update()

        self._refit_n = refit_n_temp
        self._max_nfev = max_nfev_temp
        #print(self._refit_n, self._max_nfev)
        logging.info("I found %i candidates (transitions considered: %s)." \
                     % (count, series))

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

        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        check, resol = resol_check_old(self.sess.spec, resol)
        if not check: return 0
        if self._z_off(trans_parse(series), z): return 0

        for i, s in enumerate(series.split(';')):
            self._systs_prepare()
            mod = self._syst_add(s, z, logN, b, resol)
            if mod is None: return 0
            #"""

            # Link z and b
            if i==0:
                mass_r = mass_d[trans_parse(s)[0]]
                kz = 'lines_voigt_%i_z' % mod._id
                #kb = 'lines_voigt_%i_b' % mod._id  # Not b for now
            else:
                mass = mass_d[trans_parse(s)[0]]
                self.sess.systs._constr['lines_voigt_%i_z' % mod._id] = (mod._id, 'z', kz)
                #self.sess.systs._constr['lines_voigt_%i_b' % mod._id] = \
                #    (mod._id, 'b', kb+'*%.14f' % (np.sqrt(mass_r/mass)))  # Not b for now

            if self._refit_n == 0:
                self._mods_recreate(mod_new=mod)
            self._sel_fit = True
            self._systs_cycle(mod=mod, verbose=True)
        self._spec_update()

        if recompress:
            systs._compress()

        return 0


    def systs_complete_old(self, series='all', dz=1e-4, resol=resol_def, avoid_systs=True):
        """ @brief Complete systems
        @details Add candidate transitions to fitted systems.
        @param series Series of transitions
        @param dz Threshold for redshift coincidence
        @param resol Resolution
        @param avoid_systs Avoid adding transitions over systems already fitted
        @return 0
        """
        try:
            #series = series.replace(';',',')
            #series = None if series in [None, 'None'] else str(series)
            dz = float(dz)
            resol = None if resol in [None, 'None'] else float(resol)
            avoid_systs = str(avoid_systs) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        #refit_n_temp = dc(self._refit_n)
        #max_nfev_temp = dc(self._max_nfev)

        #self._refit_n = 0
        #self._max_nfev = 1

        if series != 'all':
            t_d = trans_parse(series)
        else:
            t_d = trans_d
        #trans_arr = np.array(np.meshgrid(t_d, t_d)).T.reshape(-1,2)

        systs = dc(self.sess.systs)

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        count = 0
        #for j, syst in enum_tqdm(systs._t[3:4], len(systs._t[3:4]),
        for j, syst in enum_tqdm(systs._t, len(systs._t),
                            "cookbook_absorbers: Completing systems"):#[7:8]:
            #for i, t in enumerate(t_d):
            added = False
            for i, t in enum_tqdm(t_d, len(t_d),
                                  "cookbook_absorbers: Finding candidates"):
                s = "%s,%s" % (t,syst['series'])
                z_l, logN_l, trans_l = self.sess.lines._cands_find2(s, syst['z']-dz, syst['z']+dz, dz)
                z_len = len(z_l)-len(trans_parse(syst['series']))+1
                #print(s, z_l,trans_parse(syst['series']),trans_l)
                add = False
                if avoid_systs and 'model' in self.sess.spec.t.colnames and t in trans_l:#z_len>0:
                    for zi in z_l[0:1]:
                        for si in trans_parse(t):
                            xi = to_x(zi,si)
                            if xi > np.min(self.sess.spec.x) and xi < np.max(self.sess.spec.x):
                                wi = np.abs(self.sess.spec.x - xi).argmin()
                                mi = self.sess.spec.t['model'][wi]
                            #add = bool(add and mi>1-1e-4)
                                if mi>1-1e-4:# and False:
                                    z = z_l[0]
                                    logN = logN_l[0]
                                    add = True

                if add:
                    added = True
                    self._systs_add([t], [z], [logN], resol_list=[resol], verbose=False)
            if added:
                count += 1

        self._systs_unify()
        #self._mods_recreate()
        self._spec_update()

        if compressed:
            systs._compress()

        logging.info("I completed %i systems (transitions considered: %s)." \
                     % (count, series))
        #self._refit_n = refit_n_temp
        #self._max_nfev = max_nfev_temp


        return 0


    def systs_complete_from_z(self, series='all', series_ref=None, z_start=0,
                              z_end=6, logN=logN_def, b=b_def, resol=resol_def,
                              chi2r_thres=np.inf, dlogN_thres=np.inf,
                              refit_n=0, chi2rav_thres=1e-2,
                              max_nfev=max_nfev_def, append=True):
        """ @brief Complete systems from redshift list
        @details TBD
        @param series Series of transitions
        @param series_ref Reference series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
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
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
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

        spec = self.sess.spec
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        if series == 'all':
            trans_ex = np.unique(np.ravel([trans_parse(s)
                                           for s in systs._t['series']]))
            trans_n = list(set(trans_d)-set(trans_ex))
            #print(trans_n)
            series = ';'.join(trans_n)

        if series_ref != None:
            ws = systs._t['series']==series_ref
            wz = np.logical_and(systs._t['z']>z_start, systs._t['z']<z_end)
            w = np.logical_and(ws, wz)

            #hist, edges = np.histogram(systs._t['z'][w], bins=np.arange(0, 10, binz))
            z_list = list(systs._t['z'][w])
            id_list = list(systs._t['id'][w])
        else:
            z_list = list(systs._t['z'])
            id_list = list(systs._t['id'])


        k_list = ['lines_voigt_%i_z' % id for id in id_list]
        #print(z_list, id_list, k_list)
        for s in series.split(';'):
            s_list = [s]*len(z_list)
            logN_list = [logN]*len(z_list)
            resol_list = [resol]*len(z_list)
            self._systs_prepare(append)
            self._systs_add(s_list, z_list, logN_list, resol_list=resol_list,
                            k_list=k_list)
            self._mods_recreate()
            self._spec_update()

        if compressed:
            systs._compress()

        return 0


    def systs_improve(self, impr_n=3, refit_n=0):
        """ @brief Improve systems
        @details Improve systems adding components to reduce residuals
        @param impr_n Number of improve cycles
        @param refit_n Number of refit cycles
        @return 0
        """

        try:
            self._impr_n = int(impr_n)
            refit_n = int(refit_n)
        except:
            logging.error(msg_param_fail)
            return 0
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        logging.info("I will improve systems in at most %i iterations." \
                     % self._impr_n)
        counts = 0
        i = 0
        c = np.inf
        while i<self._impr_n and c!=0:# in range(self._impr_n):
            c = self._systs_improve()
            i += 1
            counts += c
            self.systs_fit(refit_n=refit_n)

        if compressed:
            systs._compress()

        logging.info("I improved systems in %i iterations, adding %i "
                     "components" % (i,counts))
        return 0


    def _systs_improve(self):
        """ @brief Improve systems
        @details Improve systems by adding components to reduce residuals
        @return Number of components added
        """

        spec = dc(self.sess.spec)
        lines = dc(self.sess.lines)
        systs = dc(self.sess.systs)
        ids = dc(systs._mods_t['id'])
        self.lines_find(col='deabs', append=False)

        dx_thres = np.max(spec.xmax.value-spec.xmin.value)

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        mods_sel = []
        count = 0
        #s_list = []
        #z_list = []
        for j, syst in enum_tqdm(systs._t, len(systs._t),
                            "cookbook_absorbers: Adding systems"):#[7:8]:
            l = self.sess.lines
            mod_sel = np.where([syst['id'] in i for i in ids])[0]
            #print('systs')
            #print(systs._mods_t['id'])
            #print('self.sess.systs')
            #print(self.sess.systs._mods_t['id'])
            if mod_sel not in mods_sel:
                mods_sel.append(mod_sel)
                mod = systs._mods_t['mod'][mod_sel][0]
                id = systs._mods_t['id'][mod_sel][0]
                x_mod = np.where([np.min(np.abs(x.value-mod._xf))<dx_thres
                                 for x in l.x])[0]
                if len(x_mod)>0:
                    count += 1
                    y_sel = np.argmin(l.y[x_mod].value)
                    series = np.unique([s['series'] for s in systs._t if s['id'] in id])
                    trans = np.array([])
                    for s in series:
                        trans = np.append(trans, trans_parse(s))
                    #trans = np.unique(np.ravel([trans_parse(s) for s in series]))
                    trans = np.unique(trans)
                    #print(trans)
                    z_list = [to_z(l.x[x_mod][y_sel], t) for t in trans]
                    z_sel = np.argmin([np.abs(systs._t['z']-z) for z in z_list])
                    #print([np.abs(systs._t['z']-z) for z in z_list], z_sel)
                    s = systs._t['series'][z_sel%len(systs._t)]
                    z = z_list[z_sel//len(systs._t)]
                    #s_list.append(s)
                    #z_list.append(z)
                    #print(s, z)
                    """
                    print('systs')
                    print(systs._mods_t['id'])
                    print('self.sess.systs')
                    print(self.sess.systs._mods_t['id'])
                    """
                    self._systs_add([s], [z], verbose=False)
                    """
                    print('systs')
                    print(systs._mods_t['id'])
                    print('self.sess.systs')
                    print(self.sess.systs._mods_t['id'])
                    """

            #else:
            """
            if len(m_sel)>0:
                sel = np.argmin(l.y[m_sel].value)
                s_sel = trans_parse(np.unique([systs._t['series'][np.where(systs._t['id']==i)] for i in m['id']]))
                z_sel = [np.array(systs._t['z'][systs._t['id']==i]) for i in m['id']]
                print(s_sel)
                z = [to_z(l.x[sel], s)  for s in s_sel]
                print(z)
            """
        #self._systs_add(s_list, z_list, verbose=False)
        self.sess.lines = lines
        self.lines_find(col='deabs', append=True)
        self._spec_update()

        if compressed:
            systs._compress()

        logging.info("I added %i systems." % count)
        return count


    def _abs_like(self, series='Ly-a', col='y', z_start=0, z_end=6, dz=1e-4, modul=10):
        """ @brief Assign likelihood to absorbers
        @details For each spectral bin, compute the likelihood that it has been
        absorbed by a given species.
        @param series Series of transitions
        @param col Column to apply the likelihood
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param modul Modulation of the error function
        @return likes Likelihood in redshift space for each series
        """

        try:
            #series = series.replace(';',',')
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            #resol = None if resol in [None, 'None'] else float(resol)
            modul = float(modul)
        except:
            logging.error(msg_param_fail)
            return 0

        #check, resol = resol_check_old(self.sess.spec, resol)
        #if not check: return 0

        if series == 'all':
            series = ';'.join(trans_d)

        spec = self.sess.spec
        likes = {}
        z_likes = {}

        for s in series.split(';'):
            abs = spec._t #[np.where(spec._t['y_abs'])]
            trans = trans_parse(s)
            z_all = [to_z(spec._t['x'], t) for t in trans]
            #print(z_all)
            z_int = np.arange(z_start, z_end, dz)
            #z_int = np.arange(z_end, z_start, -dz)
            z_sel, abs_sel, abs_int = [], [], []
            for z in z_all:
                sel = np.where(np.logical_and(z>z_start, z<z_end))
                if len(sel[0])>0:
                    #print(z, sel, len(sel))
                    z_sel.append(z[sel])
                    cont = spec._t['cont'][sel]
                    flux = spec._t[col][sel]
                    err = spec._t['dy'][sel]
                    #er = (cont-flux)/err/np.sqrt(2)/modul
                    #er = (np.log(flux)-np.log(cont))/(err*(-1/flux))\
                    #     /np.sqrt(2)/modul
                    #er = erf(er)

                    er = erf(np.abs((cont-flux)/err)/np.sqrt(2)/modul)

                    interp = np.interp(z_int, z[sel], er)
                    interp[z_int<np.min(z[sel])] = np.nan
                    interp[z_int>np.max(z[sel])] = np.nan
                    #print(interp)
                    abs_int.append(interp)
                    #print(np.min(z_int), np.max(z_int), np.min(z[sel]), np.max(z[sel]))

            #like = 1-np.power(1-np.nanprod(abs_int, axis=0), np.sum(~np.isnan(abs_int), axis=0))
            if len(abs_int) > 0:
                #like = 1-np.power(1-np.nanprod(abs_int, axis=0), len(trans))
                like = erfinv(np.nanprod(abs_int, axis=0))*np.sqrt(2)*modul
                like[np.isinf(like)] = -np.inf
                likes[s] = like
                z_likes[s] = z_int
                if 'like_'+s not in spec._t.colnames:
                    logging.info("I'm adding column 'like_%s' to spectrum." % s)
                    spec._t['like_'+s] = np.empty(len(spec._t))
                    spec._t['like_'+s][:] = np.nan
                else:
                    logging.warning("I'm updating column 'like_%s' in spectrum." % s)
                for t in trans:
                    """
                    if t not in spec._t.colnames:
                        logging.info("I'm adding column '%s' to spectrum." % t)
                        spec._t[t] = np.zeros(len(spec._t))
                    else:
                        logging.warning("I'm updating column '%s' in spectrum." % t)
                    """
                    if len(np.ravel([like]))==len(z_int):
                        x_int = to_x(z_int, t)
                        sel = np.logical_and(spec._t['x'].to(au.nm)>np.min(x_int),
                                             spec._t['x'].to(au.nm)<np.max(x_int))
                        cand_t = np.interp(spec._t['x'][sel].to(au.nm), x_int, like)
                        #cand_t = 1 - np.power(1-cand_t, len(trans))
                        #spec._t[t][sel] = cand_t
                        spec._t['like_'+s][sel] = np.fmax(spec._t['like_'+s][sel], cand_t)
                        #plt.step(x_int, like)

            #plt.step(z_sel[0], abs_sel[0], color='blue', alpha=0.2)
            #plt.step(z_sel[1], abs_sel[1], color='red', alpha=0.2)
            #plt.step(z_int, abs_int[0], color='blue', alpha=0.4)
            #plt.step(z_int, abs_int[-1], color='red', alpha=0.4)
            #plt.step(z_int, like)

            #plt.step(z_int, prod, color='black')
            #plt.step(spec._t['x'], cand_t, color='black')

        #plt.show()
        if 'like' not in spec._t.colnames:
            logging.info("I'm adding column 'like' to spectrum.")
            spec._t['like'] = np.empty(len(spec._t))
            spec._t['like'][:] = np.nan
        else:
            logging.warning("I'm updating column 'like' in spectrum.")
        for c in spec._t.colnames:
            if 'like_' in c:
                spec._t['like'] = np.fmax(spec._t['like'], spec._t[c])


        return likes, z_likes

    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def systs_complete_from_like(self, series='all', series_ref=None, z_start=0,
                                 z_end=6, binz=1e-2, dz=1e-4,
                                 modul=10, sigma=2, distance=3,
                                 logN=logN_def, b=b_def, resol=resol_def,
                                 chi2r_thres=np.inf, dlogN_thres=np.inf,
                                 refit_n=0, chi2rav_thres=1e-2,
                                 max_nfev=max_nfev_def, append=True):
        """ @brief Complete systems from likelihood
        @details TBD
        @param series Series of transitions
        @param series_ref Reference series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param binz Bin size to group existing redshifts
        @param dz Threshold for redshift coincidence
        @param modul Modulation of the error function
        @param sigma Threshold for accepting systems
        @param distance Distance between systems in pixels
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
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
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            binz = float(binz)
            dz = float(dz)
            modul = float(modul)
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
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


        spec = self.sess.spec
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        if series == 'all':
            trans_ex = np.unique(np.ravel([trans_parse(s)
                                           for s in systs._t['series']]))
            trans_n = list(set(trans_d)-set(trans_ex))
            #print(trans_n)
            series = ';'.join(trans_n)

        if series_ref != None:
            w = np.where(systs._t['series']==series_ref)
            #hist, edges = np.histogram(systs._t['z'][w], bins=np.arange(0, 10, binz))

            # Rebecca's fix
            median_z = np.nanmedian(systs._t['z'][w])
            z_lower_boundary = np.round(((median_z + 0.5 * binz) % binz) / dz) * dz
            hist, edges = np.histogram(systs._t['z'][w], bins=np.arange(z_lower_boundary, 10, binz))
        else:
            hist, edges = np.histogram(systs._t['z'], bins=np.arange(0, 10, binz))

        w_z = np.logical_and(edges>z_start, edges<z_end)


        self._likes, self._z_likes = {}, {}
        for z in edges[:-1][hist>0]:
            if z>z_start and z+binz<z_end:
                z_s = z
                z_e = z+binz
            elif z<z_start and z+binz>z_end:
                z_s = z_start
                z_e = z_end
            elif z<z_start and z+binz>z_start:
                z_s = z_start
                z_e = z+binz
            elif z<z_end and z+binz>z_end:
                z_s = z
                z_e = z_end
            else:
                z_s = np.nan
                z_e = np.nan
            #print(z, z+binz, z_start, z_end, z_s, z_e)
            if not np.isnan(z_s) and not np.isnan(z_e):
                likes, z_likes = self._abs_like(series, 'y', z_s, z_e, dz, modul)
                """
            for s in likes.keys():
                if s not in self._likes.keys():
                    self._likes[s] = likes[s]
                    self._z_likes[s] = z_likes[s]
                else:
                    self._likes[s] = np.append(self._likes[s], likes[s])
                    self._z_likes[s] = np.append(self._z_likes[s], likes[s])
                """
                self._likes = likes
                self._z_likes = z_likes
                self._systs_like(series, sigma, distance, logN, b, resol, chi2r_thres,
                                 dlogN_thres, refit_n, chi2rav_thres, max_nfev, append)

        if compressed:
            systs._compress()
        """
        logging.info("I'm completing systems with %i additional transitions." \
                     % len(self._likes))
        self._systs_like(series, thres, logN, b, resol, chi2r_thres,
                         dlogN_thres, refit_n, chi2rav_thres, max_nfev, append)
        """
        return 0


    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def systs_new_from_like(self, series='Ly-a', col='y', z_start=0, z_end=6,
                            dz=1e-4, modul=10, sigma=2, distance=3,
                            logN=logN_def, b=b_def, resol=resol_def,
                            chi2r_thres=np.inf, dlogN_thres=np.inf,
                            refit_n=0, chi2rav_thres=1e-2, max_nfev=max_nfev_def,
                            append=True):

        """ @brief New systems from likelihood
        @details Add Voigt models by testing transitions in a given redshift
        range and assignign the most likely identification to absorption
        features.
        @param series Series of transitions
        @param col Column to apply the likelihood
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param modul Modulation of the error function
        @param sigma Threshold for accepting systems
        @param distance Distance between systems in pixels
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
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
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            modul = float(modul)
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
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

        check, resol = resol_check_old(self.sess.spec, resol)
        if not check: return 0

        self._likes, self._z_likes = self._abs_like(series, col, z_start, z_end, dz,
                                                    modul)
        self._systs_like(series, sigma, distance, logN, b, resol, chi2r_thres,
                         dlogN_thres, refit_n, chi2rav_thres, max_nfev, append)

        return 0


    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def _systs_like(self, series='Ly-a', sigma=2, distance=3, logN=logN_def,
                    b=b_def, resol=resol_def, chi2r_thres=np.inf,
                    dlogN_thres=np.inf, refit_n=0, chi2rav_thres=1e-2,
                    max_nfev=max_nfev_def, append=True):

        try:
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
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

        if series == 'all':
            series = ';'.join(trans_d)

        spec = self.sess.spec
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()


        likes = self._likes
        z_likes = self._z_likes
        series_split = series.split(';')

        def _nan_check(spec, s_list, z_list, logN_list, resol_list):
            # Check that components do not fall in masked regions
            sel = []
            for ssel,zsel in zip(s_list, z_list):
                ysel = []
                for t in trans_parse(ssel):
                    xsel = to_x(zsel, t)
                    ysel.append(np.interp(xsel, spec._t['x'],
                                          spec._t['y']))
                sel.append(not np.any(np.isnan(ysel)))
            wsel = np.where(sel)[0]
            s_list = np.array(s_list)[wsel]
            z_list = np.array(z_list)[wsel]
            logN_list = np.array(logN_list)[wsel]
            resol_list = np.array(resol_list)[wsel]
            return s_list, z_list, logN_list, resol_list



        k_list = []
        id_list = []
        for i, s in enumerate(series_split):
            #print(s, 'systs_like')
            #series_o = list(set(series_split) - set([s]))
            trans = trans_parse(s)
            #z_int = np.arange(z_start, z_end, dz)
            if s in likes.keys():
                #print(likes[s])
                z_int = z_likes[s]
                #print(s)
                #plt.plot(z_int, likes[s])
                #w = np.where(likes[s]>sigma)
                #plt.plot(z_int[w], likes[s][w])
                #print(len(w[0]))
                """
                for s_o in series_o:
                    trans_o = trans_parse(s_o)
                    for t in trans_o:
                        print(t)
                        t_o = np.interp(x_w, spec._t['x'], spec._t[t])
                        print(likes_o)
                """
                #p0, _ = find_peaks(likes[s][w], distance=distance)
                #plt.scatter(z_int[w][p0], likes[s][w][p0])

                p0, _ = find_peaks(likes[s], distance=distance)
                pw = np.where(likes[s][p0]>sigma)
                p0 = p0[pw]
                #plt.scatter(z_int[p0], likes[s][p0])

                """
                # Check if likelihood peaks are higher than those of all other
                # transitions at those wavelengths
                x_w = np.array([to_x(z_int[w][p0], t) for t in trans])
                #for x_wi in x_w:
                #    plt.scatter(x_wi, likes[s][w][p0])
                t_all = np.array([])
                #for so in np.array(series_split):
                for so in likes.keys():
                    t_all = np.append(t_all, trans_parse(so))
                #t_all = np.ravel([trans_parse(so) for so in np.array(series_split)])
                liket = np.array([])
                for to in t_all:
                    liket = np.append(liket, np.interp(x_w, spec._t['x'], spec._t[to]))
                    #print(to, liket)
                liket = np.reshape(liket, (len(t_all), len(np.ravel(x_w))))
                liket[np.isnan(liket)] = -9999
                check = [t in trans for t in t_all[np.argmax(liket, axis=0)]]
                sel = np.prod(np.reshape(check, (len(trans), len(p0))), axis=0)
                #print(sel)
                #print(p0)
                #print(p0[np.where(sel)])
                """
                #print(z_int[w][p0])
                #x_w = to_x(z_int[w][p0], trans[0])
                x_w = to_x(z_int[p0], trans[0])
                #print(x_w)
                psel = []
                for p, x in zip(p0, x_w):
                    c = np.abs(x - spec._t['x']).argmin()
                    if spec._t['like_'+s][c]==spec._t['like'][c]:
                        psel.append(p)


                p = psel

                for ssub in s.split(':'):
                    if k_list == []:
                        s_list = [ssub]*len(p)

                        #z_list = z_int[w][p]
                        z_list = z_int[p]
                        logN_list = [logN]*len(p)
                        resol_list = [resol]*len(p)

                        # Check that components do not fall in masked regions
                        """
                        sel = []
                        for ssel,zsel in zip(s_list, z_list):
                            xint = []
                            for t in trans_parse(ssel):
                                xsel = to_x(zsel, t)
                                xint.append(np.interp(xsel, spec._t['x'],
                                                      spec._t['y']))
                            sel.append(not np.any(np.isnan(xint)))
                        wsel = np.where(sel)[0]
                        s_list = np.array(s_list)[wsel]
                        z_list = np.array(z_list)[wsel]
                        logN_list = np.array(logN_list)[wsel]
                        resol_list = np.array(resol_list)[wsel]
                        """
                        s_list, z_list, logN_list, resol_list = \
                            _nan_check(spec, s_list, z_list, logN_list, resol_list)

                        if len(s_list)>0:
                            self._systs_prepare(append)
                            id_list = self._systs_add(s_list, z_list, logN_list,
                                                      resol_list=resol_list)
                            self._spec_update()
                        else:
                            id_list = []
                    else:
                        s_list = [ssub]*len(z_list)
                        logN_list = [logN]*len(z_list)
                        resol_list = [resol]*len(z_list)

                        s_list, z_list, logN_list, resol_list = \
                            _nan_check(spec, s_list, z_list, logN_list, resol_list)

                        self._systs_prepare(append)
                        id_list = self._systs_add(s_list, z_list, logN_list,
                                                  resol_list=resol_list, k_list=k_list)
                        self._mods_recreate()
                        self._spec_update()
            else:
                id_list = []
            if i == 0:
                k_list = ['lines_voigt_%i_z' % id for id in id_list]
                #print(k_list)
        #plt.show()
        if compressed:
            systs._compress()
        #print('end')
        return 0


    def systs_merge(self, to_row=0, from_rows=[1]):
        """ @brief Merge a system into the current system
        @details Merged systems appear as a single entry in the compressed
        system table.
        @param num1 Row of the current system
        @param num2 Row of the system to be merged
        @return 0
        """

        try:
            to_row = int(to_row)-1
            from_rows = np.array(from_rows)-1
        except:
            logging.error(msg_param_fail)
            return 0

        t = self.sess.systs.t
        #print(t['z'][from_rows])
        z_app = np.append(t['z'][from_rows], t['z'][to_row])
        logN_app = np.append(t['logN'][from_rows], t['logN'][to_row])

        t['z'][to_row] = np.average(z_app, weights=10**logN_app)
        t['logN'][to_row] = np.log10(np.sum(10**logN_app))
        t.remove_rows(from_rows)
        t.sort(['z','id'])

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
        @param b Guess Doppler broadening
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

        check, resol = resol_check_old(self.sess.spec, resol)
        if not check: return 0

        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()


        for s in series.split(';'):
            z_list, y_list = self._lines_cands_find(s, z_start, z_end, dz)
            z_list, logN_list, _ = self.sess.lines._cands_find2(s, z_start, z_end,
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

        if compressed:
            systs._compress()

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
        @param b Guess Doppler broadening
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

    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def _systs_new_from_erf(self, series='Ly-a', col='y', z_start=0, z_end=6,
                            sigma=1, distance=3, append=True):
        """ @brief New systems from error function
        @details TBD
        @param series Series of transitions
        @param col Column to apply the likelihood
        @param z_start Start redshift
        @param z_end End redshift
        @param sigma Significance of absorbers (in units of the local error)
        @param distance Distance between systems in pixels
        @param append Append systems to existing system list
        @return 0
        """

        try:
            z_start = float(z_start)
            z_end = float(z_end)
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
        except:
            logging.error(msg_param_fail)
            return 0

        modul = 10
        distance = 3
        self.systs_new_from_like(series=series, col=col, z_start=z_start,
                                 z_end=z_end, modul=modul, sigma=sigma,
                                 distance=distance, append=append)

        return 0

    def _series_fit(self, series, zem, z_start=None, z_end=None, sigma=2, iter_n=3):

        def z_check(zem, z_start, z_end, s):
            if zem != None:
                if 'Ly_a' not in trans_parse(s):
                    z_start = (1+zem)*xem_d['Ly_a']/xem_d[trans_parse(s)[-1]]-1
                else:
                    z_start = (1+zem)*xem_d['Ly_b']/xem_d['Ly_a']-1
                z_end = zem
            return z_start, z_end

        def resid_peaks(std=10):
            spec = self.sess.spec
            #spec._gauss_convolve(std=std, input_col='deabs', output_col='deabs_mod')
            spec.t['deabs_mod'] = np.maximum(1, spec.t['deabs'])
            p, _ = find_peaks(spec.t['deabs_mod'], prominence=0.05)
            #print(np.array(spec._t['x'][p]))

        for s in series.split(';'):
            z_start, z_end = z_check(zem, z_start, z_end, s)
            self._systs_new_from_erf(series=s, z_start=z_start, z_end=z_end,
                                     sigma=sigma)
        self.systs_fit(refit_n=1)#, max_nfev=0)
        for i in range(iter_n):
            self.sess.systs._freeze_pars()
            for s in series.split(';'):
                z_start, z_end = z_check(zem, z_start, z_end, s)
                self._systs_new_from_erf(series=s, col='deabs',
                                         z_start=z_start, z_end=z_end, sigma=sigma)
            self.systs_fit(refit_n=0)#, max_nfev=0)
            self.sess.systs._unfreeze_pars()
            resid_peaks()
        if iter_n > 0:
            self.systs_fit(refit_n=1)
        #plt.show()
        return 0


    def lya_fit(self, zem=None, z_start=None, z_end=None, sigma=1, iter_n=3):
        """ @brief Fit the Lyman-alpha forest
        @details The recipe identifies Lyman-alpha absorbers using the
        likelihood method and fits them. The procedure is iterated on residuals
        of the fit to improve it.
        @param zem Emission redshift
        @param z_start Start redshift (ignored if zem is specified)
        @param z_end End redshift (ignored if zem is specified)
        @param sigma Significance of absorbers (in units of the local error)
        @param iter_n Number of iterations on residuals
        @return 0
        """

        try:
            zem = None if zem in [None, 'None'] else float(zem)
            z_start = None if z_start in [None, 'None'] else float(z_start)
            z_end = None if z_end in [None, 'None'] else float(z_end)
            sigma = float(sigma)
            iter_n = int(iter_n)
        except:
            logging.error(msg_param_fail)
            return 0

        if zem != None:
            z_start = (1+zem)*xem_d['Ly_b']/xem_d['Ly_a']-1
            z_end = zem

        self._series_fit('Ly_a', zem, z_start, z_end, sigma, iter_n)
        #plt.show()
        return 0


    def lyab_chunk_fit(self, z_start, z_end, dv=1000, ol=100):
        """ @brief Fit the Lyman-alpha and Lyman-beta forest by chunks
        @details The recipe fits previously detected Lyman-alpha and Lyman-beta
        absorbers by chunks. The chunks span from `z_start` to `z_end` with a
        width `dv` and an overlap `ol`.
        @param z_start Start redshift
        @param z_end End redshift
        @param dv Width of the chunk (km/s)
        @param ol Overlap of adjacent chunks (km/s)
        @return 0
        """

        try:
            z_start = float(z_start)
            z_end = float(z_end)
            dv = float(dv)
            ol = float(ol)
        except:
            logging.error(msg_param_fail)
            return 0

        z_range = []
        z = z_start
        while z<z_end:
            dz = dv/aconst.c.to(au.km/au.s).value*(1+z)
            dzo = (dv-ol)/aconst.c.to(au.km/au.s).value*(1+z)
            z_range.append((z,z+dz))
            z += dzo

        chunks = [(i, i+1) for i in range(571,583)]
        for z_c in z_range:
            c_lyb = (1+np.array(z_c)) * xem_d['Ly_b'].value
            c_lya = (1+np.array(z_c)) * xem_d['Ly_a'].value
            chunks = "%s-%s;%s-%s" % (c_lyb[0], c_lyb[1], c_lya[0], c_lya[1])
            self.chunk_fit(chunks, recreate=False)
        self._mods_recreate()

        return 0


    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def red_fit(self,
                series='CIV;SiIV:CIV;SiII_1526:CIV;AlIII;MgII_2796,MgII_2803;'\
                       +'FeII_2382,FeII_2600:MgII_2796,MgII_2803',
                zem=None, z_start=None, z_end=None, sigma=1, iter_n=3):
        """ @brief Fit the red part of the spectrum forest
        @details The recipe identifies Lyman-alpha absorbers using the
        likelihood method and fits them. The procedure is iterated on residuals
        of the fit to improve it.
        @param zem Emission redshift
        @param z_start Start redshift (ignored if zem is specified)
        @param z_end End redshift (ignored if zem is specified)
        @param sigma Significance of absorbers (in units of the local error)
        @param iter_n Number of iterations on residuals
        @return 0
        """

        try:
            zem = None if zem in [None, 'None'] else float(zem)
            z_start = None if z_start in [None, 'None'] else float(z_start)
            z_end = None if z_end in [None, 'None'] else float(z_end)
            sigma = float(sigma)
            iter_n = int(iter_n)
        except:
            logging.error(msg_param_fail)
            return 0

        """
        for s in ['CIV','MgII_2796,MgII_2803','SiIV']:
            t = trans_parse(s)
            if zem != None:
                z_start = (1+zem)*xem_d['Ly_a']/xem_d[t[-1]]-1
                z_end = zem
        """
        self._series_fit(series, zem, z_start, z_end, sigma, iter_n)
            #self._systs_new_from_erf(series=s, z_start=z_start, z_end=z_end,
            #                         sigma=sigma)

        return 0


    def systs_complete(self, x, x_doublet="", dz_systs=1e-4, dz_doublet=1e-5,
                       dz_unknown=1e-5, dz_galact=1e-2, z_min=0, z_max=9):
        """ @brief Complete system identification
        @details Complete system identification by (1) associating unknown lines
        to known systems; (2) identifying doublets from their wavelength ratio;
        (3) identifying unknown lines by finding possible redshift coincidences;
        and (4) identifying possible galactic lines. Identification is done
        using a reference list of transitions)
        @param x List of line wavelengths (nm)
        @param x_doublet Test wavelengths of the doublets (nm; e.g. 500,501;600-601)
        @param dz_systs Threshold for redshift coincidence with known systems
        @param dz_doublet Threshold for redshift coincidence between doublet members
        @param dz_unknown Threshold for redshift coincidence between unknown lines
        @param dz_galact Threshold for coincidence with redshift 0
        @param z_min Minimum redshift
        @param z_max Maximum redshift
        @return 0
        """

        try:
            x = np.array(x[1:-1].split(','), dtype=float) if type(x)==str \
                else np.array(x)
            dz_systs = float(dz_systs)
            dz_doublet = float(dz_doublet)
            dz_unknown = float(dz_unknown)
            dz_galact = float(dz_galact)
            z_min = float(z_min)
            z_max = float(z_max)
        except:
            logging.error(msg_param_fail)
            return 0

        self._systs_assoc(x, dz_systs)
        if x_doublet not in [None, ""]: self._doublet_iden(x_doublet, dz_doublet)
        self._unknown_iden(x, dz_unknown, z_min, z_max)
        self._galact_iden(x, dz_galact)

        return 0


    def _systs_assoc(self, x, dz=1e-4):
        """ @brief Associate lines to systems
        @details Associate unknown lines to known systems by finding possible
        redshift coincidences (using a reference list of transitions)
        @param x List of line wavelengths (nm)
        @param dz Threshold for redshift coincidence
        @return 0
        """

        try:
            x = np.array(x[1:-1].split(','), dtype=float) if type(x)==str \
                else np.array(x)
            eps = float(dz)
        except:
            logging.error(msg_param_fail)
            return 0

        #print('\n\n-----------------------------------------------------\nATOM_PAR TABLE\n-----------------------------------------------------\n',atom_par[0:226])
        wl = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col2']]))
        labels = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col1']]))
        t = self.sess.systs._t
        #print('\n\n-----------------------------------------------------\nSYSTEM TABLE\n-----------------------------------------------------\n',t)
        z_systs = np.ravel(np.array([d for d in t['z']]))
        corr = np.ravel(np.array([d for d in t['series']]))

        x.sort()

        #print('\n\n------------------------------------------------\nMATCH UNKOWN LINES TO KNOWN SYSTEMS\n------------------------------------------------\n')
        for l in x:
            #print(f'\n------------- line @ {l} -------------')
            #print('Ion\t\tref\t\tz')

            z = l/wl-1
            #eps = 1e-4

            lab_l = []
            corr_l = []
            z_ref_l = []
            for i, lab in enumerate(labels):
                for j,z_ref in enumerate(z_systs):
                    dz = z[i]-z_ref

                    if (np.abs(dz)<eps) and (corr[j]!='Ly_a'):
                        lab_l.append(lab)
                        corr_l.append(corr[j])
                        z_ref_l.append(z_ref)
                        #print(f'{lab}\t{corr[j]}\t{z_ref}')
            lab_lu = np.unique(lab_l)
            if len(lab_lu)>0:
                logging.info("Line at %3.4f has %i possible association%s:" \
                            % (l, len(lab_lu), 's' if len(lab_lu)>1 else ''))
                #print(lab_lu, lab_l, corr_l, z_ref_l)
                for u in lab_lu:
                    w = np.where(np.array(lab_l)==u)[0]
                    logging.info(" %s at z=%3.4f (as %s)" \
                                 % (', '.join(np.array(corr_l)[w]), np.mean(np.array(z_ref_l)[w]), u))


    def _doublet_iden(self, x, dz=1e-5):
        """ @brief Identify doublets
        @details Identify doublets from their wavelength ratio (using a
        reference list of transitions)
        @param x Test wavelengths of the doublets (nm; e.g. 500,501;600-601)
        @param dz Threshold for redshift coincidence
        @return 0
        """

        try:
            x = doublet_parse(x)
            eps = float
        except:
            logging.error(msg_param_fail)
            return 0

        #print('\n\n-----------------------------------------------------\nATOM_PAR TABLE\n-----------------------------------------------------\n',atom_par[0:226])
        wl = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col2']]))
        labels = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col1']]))
        t = self.sess.systs._t

        x.sort()

        for (l1, l2) in x:
            if l1 > l2:
                l1, l2 = l2, l1
            #print('\n\n------------------------------------------------\nDO TWO UNIDENTIFIED LINES HAVE A SIMILAR SHAPE?\n------------------------------------------------\n')
            ratio = l1/l2
            eps = 1e-5
            for i in range(len(wl)):
                for j in range(i, len(wl)):
                    ratio_i = wl[i]/wl[j]

                    z = l1/wl[i]-1
                    if (np.abs(ratio_i-ratio)<=eps):
                        #print('----------------------------')
                        #print(labels[i],'\t',wl[i])
                        #print(labels[j],'\t',wl[j])
                        #print(f'proposed z = {l1/wl[i]-1}')
                        logging.info("Doublet at %3.4f-%3.4f nm is compatible "\
                                     "with %s-%s at z=%3.4f."
                                     % (l1,l2,labels[i],labels[j],l1/wl[i]-1))


    def _unknown_iden(self, x, dz=1e-5, z_min=0, z_max=9):
        """ @brief Identify unknown lines
        @details Identify unknown lines by finding possible redshift
        coincidences between them (using a reference list of transitions)
        @param x List of line wavelengths (nm)
        @param dz Threshold for redshift coincidence
        @param z_min Minimum redshift
        @param z_max Maximum redshift
        @return 0
        """

        try:
            x = np.array(x[1:-1].split(','), dtype=float) if type(x)==str \
                else np.array(x)
            eps = float(dz)
            z_min = float(z_min)
            z_max = float(z_max)
        except:
            logging.error(msg_param_fail)
            return 0

        #print('\n\n-----------------------------------------------------\nATOM_PAR TABLE\n-----------------------------------------------------\n',atom_par[0:226])
        wl = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col2']]))
        labels = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col1']]))
        t = self.sess.systs._t

        x.sort()

        #print('\n\n------------------------------------------------\nWAVELENGTH RATIOS BETWEEN ALL THE UNKNOWN LINES\n------------------------------------------------\n')
        #eps = 1e-5
        #z_min, z_max = -1e-4, 4

        skip = np.array(['H','Ly','D'])

        iden_l = []
        lab_l = []
        z_l = []
        for i,_ in enum_tqdm(range(len(x)), len(x),
                             "cookbook_absorbers: Identifying unknown lines"):
            for j in range(i+1,len(x)):
                ratio = x[i]/x[j]
                for k in range(len(wl)):
                    for l in range(k+1, len(wl)):
                        ratio_em = wl[k]/wl[l]
                        z1 = x[i]/wl[k]-1
                        z2 = x[j]/wl[l]-1
                        if ((np.abs(ratio_em-ratio)<=eps) and (np.abs(z1-z2)<=eps) and (z1<=z_max) and (z1>=z_min)):
                            if np.array([N not in labels[k] for N in skip]).all and np.array([N not in labels[l] for N in skip]).all():
                                #print('----------------------------------------')
                                #print(f'line 1: {x[i]:.2f}, {labels[k]}')
                                #print(f'line 2: {x[j]:.2f}, {labels[l]}')
                                #print(f'proposed z = {z1:.4f}')
                                iden_l.append(x[i])
                                iden_l.append(x[j])
                                lab_l.append(labels[k])
                                lab_l.append(labels[l])
                                z_l.append(z1)
                                z_l.append(z1)

        iden_lu = np.unique(iden_l)
        for i in iden_lu:
            w = np.where(np.array(iden_l)==i)[0]
            logging.info("Line at %3.4f has %i possible redshift coincidence%s:" \
                        % (i, len(w), 's' if len(w)>1 else ''))
            z = np.array(z_l)[w]
            lab = np.array(lab_l)[w]
            zs = z[np.argsort(z)]
            labs = lab[np.argsort(z)]
            for (zsi, labsi) in zip(zs, labs):
                wz = np.where(np.logical_and(np.array(z_l)==zsi,
                                             np.array(lab_l)!=labsi))[0]
                #print(wz)
                #print(iden_lu)
                logging.info(" z=%3.4f (as %s, with %s at %3.4f)"  \
                             % (zsi, labsi, np.array(lab_l)[wz][0],
                                np.array(iden_l)[wz][0]))

        #print('\n\n--------------------------------------------\nANY GALACTIC TRANSITIONS?\n--------------------------------------------\n')

    def _galact_iden(self, x, dz=1e-2):
        """ @brief Identify galactic lines
        @details Identify unknown lines by finding possible coincidences at
        redshift 0 with a reference list of transitions
        @param x List of line wavelengths (nm)
        @param dz Threshold for redshift coincidence
        @return 0
        """

        try:
            x = np.array(x[1:-1].split(','), dtype=float) if type(x)==str \
                else np.array(x)
            eps = float(dz)
        except:
            logging.error(msg_param_fail)
            return 0

        #print('\n\n-----------------------------------------------------\nATOM_PAR TABLE\n-----------------------------------------------------\n',atom_par[0:226])
        wl = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col2']]))
        labels = np.ravel(np.array([d2 for d2 in atom_par[0:226]['col1']]))
        t = self.sess.systs._t

        x.sort()

        l_l = []
        lab_l = []
        z_l = []
        for i, l in enumerate(x):
            #print(f'\n------------- line @ {l} -------------')
            #print('Ion\t\tz')

            z = l/wl-1
            #eps = 1e-2
            for j, lab in enumerate(labels):
                if abs(z[j])<eps:
                    #print(f'{lab}\t{z[j]:.4f}')
                    l_l.append(l)
                    lab_l.append(lab)
                    z_l.append(z[j])

        l_lu = np.unique(l_l)
        for l in l_lu:
            w = np.where(np.array(l_l)==l)[0]
            logging.info("Line at %3.4f has %i possible galactic identification%s:" \
                        % (i, len(w), 's' if len(w)>1 else ''))
            lab = np.array(lab_l)[w]
            z = np.array(z_l)[w]
            zs = z[np.argsort(z)]
            labs = lab[np.argsort(z)]
            for (zsi, labsi) in zip(zs, labs):
                logging.info(" z=%3.4f (as %s)" % (zsi, labsi))