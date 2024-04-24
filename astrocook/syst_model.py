from .functions import *
from .spectrum import Spectrum
from .vars import *
from astropy import table as at
from copy import deepcopy as dc
import datetime
from lmfit import CompositeModel as LMComposite
from lmfit import Model as LMModel
from lmfit import Parameters as LMParameters
from matplotlib import pyplot as plt
import numpy as np
import operator
import time
import warnings
warnings.filterwarnings("ignore")

# from line_profiler import LineProfiler
# from filprofiler.api import profile

class SystModel(LMComposite):

    def __init__(self, systs, series=[], vars=None, constr=None, z0=None,
                 lines_func=lines_voigt,
                 psf_func=None,
                 cont_func=None):
        #self._spec = spec
        try:
            self._mods_t = systs._mods_t
        except:
            self._mods_t = None
        self._id = systs._id
        self._series = series
        if vars is None: vars = {}
        if constr is None: constr = {}
        self._vars = vars
        self._constr = constr
        self._z0 = z0
        self._lines_func = lines_func
        self._psf_func = psf_func
        self._active = True

    def _fit(self, fit_kws={}):
        vary = np.any([self._pars[p].vary for p in self._pars])
        fit_kws_c = dc(fit_kws)
        if vary:
            time_start = datetime.datetime.now()
            if 'use_jac' in fit_kws_c:
                use_jac = fit_kws_c['use_jac']
                del fit_kws_c['use_jac']
            else:
                use_jac = False
            if 'max_nfev' in fit_kws_c:
                max_nfev = fit_kws_c['max_nfev']
                del fit_kws_c['max_nfev']
            else:
                max_nfev = None

            plot_jac = False
            if plot_jac: col = 1

            def _jac(x0, *args, **kwargs):
                return globals()[self._lines_func.__name__+'_jac']\
                    (x0, self._xf, series=self._series, resol=self._resol,
                    spec=self._spec, *args, **kwargs)
            pars = [self._pars[p].value for p in self._pars if 'z' in p
                    or 'logN' in p or 'b' in p and 'btur' not in p]
            if plot_jac: plt.plot(self._xf, _jac(pars)[:,col], color='red')

            plt.plot(self._xf, self._yf)
            if use_jac:
                fit_kws_c['jac'] = _jac

                try:
                    fit = super(SystModel, self).fit(self._yf, self._pars, x=self._xf,
                                                     weights=self._wf,
                                                     max_nfev=max_nfev,
                                                     fit_kws=fit_kws_c,
                                                     nan_policy='omit',
                                                     method='least_squares')
                    if plot_jac: plt.plot(self._xf, fit.jac[:,col], color='green')
                except:
                    del fit_kws_c['jac']
                    use_jac = False
            if not use_jac:
                fit = super(SystModel, self).fit(self._yf, self._pars, x=self._xf,
                                                 weights=self._wf,
                                                 max_nfev=max_nfev,
                                                 fit_kws=fit_kws_c,
                                                 nan_policy='omit',
                                                 method='least_squares')
                #class_unmute(self, Spectrum)
                if plot_jac: plt.plot(self._xf, fit.jac[:,col], color='blue')
            if plot_jac: plt.show()
            time_end = datetime.datetime.now()
            self._pars = fit.params
            self._ys = self.eval(x=self._xs, params=self._pars)
            #plt.plot(self._xs, self._ys)
            #plt.plot(self._xs[prova], self.eval(x=self._xs[prova], params=self._pars))
            #plt.plot(self._xf, self.eval(x=self._xf, params=self._pars))
            #plt.show()
            self._chi2r = fit.redchi
            self._aic = fit.aic
            self._bic = fit.bic
            return 0
        else:
            return 1


    def _make_comp(self):
        super(SystModel, self).__init__(self._group, self._psf, convolve_simple)
        #self._pars.pretty_print()

    """
    def _zero(self, x):
        return 0*x
    """

    def _make_comp2(self):
#        super(SystModel, self).__init__(self._group, LMModel(self._zero), operator.add)
#        super(SystModel, self).__init__(self._group, LMModel(zero), operator.add)
        if hasattr(self, '_group'):
            super(SystModel, self).__init__(self._group, LMModel(zero), operator.add)
        #self._pars.pretty_print()


    def _make_defs(self, defs=None):
        if defs is None:
            self._defs = dc(pars_std_d)
        else:
            defs_complete = {}
            for d in pars_std_d:
                if d in defs:
                    defs_complete[d] = defs[d]
                else:
                    defs_complete[d] = pars_std_d[d]
            self._defs = defs_complete
        for v in self._vars:
            if v in self._defs:
                self._defs[v] = self._vars[v]
        if 'x_lim' in defs:
            if defs['x_lim'] is None:
                self._defs['x_lim'] = None
            else:
                xl = defs['x_lim']
                self._defs['x_lim'] = [[float(i) for i in x.split('-')] for x in xl.split(',')]


    def _make_groups(self, spec, thres):
        """ @brief Group lines that must be fitted together into a single model.
        """

        time_check = False
        #T = time.time()
        #tt = time.time()
        #spec = self._spec

        mods_t = self._mods_t
        d = self._defs

        self._xs = np.array(spec._safe(spec.x).to(au.nm))
        if d['x_lim'] is None:
            self._xl = np.array(spec._safe(spec.x).to(au.nm))
        else:
            x = spec._safe(spec.x).to(au.nm).value
            w = np.array([], dtype=int)
            for xl in d['x_lim']:
                w = np.append(w, np.where(np.logical_and(x>xl[0], x<xl[1]))[0])
            w = np.unique(w)
            self._xl = x[w]

        ys = self._lines.eval(x=self._xs, params=self._pars)
        yl = self._lines.eval(x=self._xl, params=self._pars)
        if np.abs(np.min(ys)-np.min(yl))>thres and np.min(yl) > 1-thres:
            self._active = False

        ysmin = np.amin(ys)
        self._group = self._lines
        self._group_list = []
        modified = False
        #print(mods_t['id'])
        #if time.time()-tt > 1e-2: print('A %.4f' % (time.time()-tt), end=' ')
        #tt = time.time()
        for i, s in enumerate(mods_t):
            modified = False

            if time_check:
                ttt = time.time()
            #print(s['id'])
            mod = s['mod']
            ys_s = mod._ys

            if time_check:
                print('0 %.4f' % (time.time()-ttt))
                ttt = time.time()
            cond, pars_cond = False, False
            for p,v in self._constr.items():
                for mod_p,mod_v in mod._pars.items():
                    cond = cond or v.split('*')[0]==mod_p
            if cond: pars_cond = True
            if time_check:
                print('1 %.4f' % (time.time()-ttt))
                ttt = time.time()
            if not cond:
            #print(s['id'])
                #ymax = np.maximum(ys, ys_s)
                yminmax = np.amin(np.maximum(ys,ys_s))
                ymin = np.amin(ys)
                if time_check:
                    print('2 %.4f' % (time.time()-ttt))
                    ttt = time.time()
                #cond = np.amin(ymax)<1-thres or np.amin(ymax)==np.amin(ys)
                cond = yminmax<1-thres*(1-ymin) or yminmax==ysmin
                if time_check:
                    print('3 %.4f' % (time.time()-ttt))
                    ttt = time.time()

            if time_check:
                print('4 %.4f' % (time.time()-ttt))
                ttt = time.time()
            if cond and self._active==mod._active: #y_cond or pars_cond:
                #print('mod')
                #mod._pars.pretty_print()
                #print('self')
                #self._pars.pretty_print()
                #try:
                if time_check:
                    print('b %.4f' % (time.time()-ttt))
                    ttt = time.time()

                self._group *= mod._group
                if time_check:
                    print('a %.4f' % (time.time()-ttt))
                    ttt = time.time()
                #except:
                #    self._group_sel = -1
                #    return
                mod = s['mod']
                for p,v in mod._pars.items():
                    if v.expr != None:
                        self._constr[p] = v.expr
                        v.expr = ''
                self._pars.update(mod._pars)
                if time_check:
                    print('5 %.4f' % (time.time()-ttt))
                    ttt = time.time()
                if pars_cond or self._constr != {}:
                    for p,v in self._constr.items():
                        if all([p in self._pars for p in v.split('+10**')[1:]]):
                            if 'backup' not in p:
                                self._pars[p].set(expr = v)
                        if v.split('+')[0] in self._pars:
                            if 'backup' not in p:
                                self._pars[p].set(expr = v)
                        #"""
                        #self._pars[p].expr = v
                        if v != '' and '**' not in v:
                            vs = v.split('*')
                            f = float(vs[1]) if len(vs)==2 else 1
                            try:
                                self._pars[p].min = self._pars[vs[0]].min
                                self._pars[p].max = self._pars[vs[0]].max
                                self._pars[p].value = self._pars[vs[0]].value * f
                            except:
                                pass
                        """
                                vs = v.split('+')
                                f = float(vs[1]) if len(vs)==2 else 1
                                try:
                                    self._pars[p].min = self._pars[vs[0]].min
                                    self._pars[p].max = self._pars[vs[0]].max
                                    self._pars[p].value = self._pars[vs[0]].value + f
                                except:
                                    self._pars[p].expr = ''
                        """
                self._group_list.append(i)
                if time_check:
                    print('6 %.4f' % (time.time()-ttt))
                    ttt = time.time()
                #mod_ys = self._group.eval(x=self._xs, params=self._pars)
                mod._ys = ys_s*ys
                ys = mod._ys
                #print('mid', np.mean(np.abs(mod_ys-mod._ys)))
                modified = True
                if time_check:
                    print('7 %.4f' % (time.time()-ttt))
                    ttt = time.time()
                #print('')
                #print(i, id(self._group), id(self._pars))
            if time_check:
                print('8 %.4f' % (time.time()-ttt))

        #if time.time()-tt > 1e-2: print('B %.4f' % (time.time()-tt), end=' ')
        #tt = time.time()
        if len(self._group_list) > 1:
            ids = [i for il in mods_t['id'][self._group_list[1:]] for i in il]
            mods_t.remove_rows(self._group_list[1:])
            for i in ids:
                mods_t['id'][self._group_list[0]].append(i)
        #if time.time()-tt > 1e-2: print('b %.4f' % (time.time()-tt), end=' ')
        #tt = time.time()
        if self._group_list == []:
            self._group_sel = -1
        else:
            self._group_sel = self._group_list[0]
        #self._ys = self._group.eval(x=self._xs, params=self._pars)
        if modified: # and False:
            self._ys = mod._ys
            #print('hey')
        else:
            #self._ys = self._group.eval(x=self._xs, params=self._pars)
            self._ys = ys
            """
            try:
                print('end', modified, np.mean(np.abs(self._ys-ys)))
            except:
                pass
            """
        #if time.time()-tt > 1e-2: print('C %.4f' % (time.time()-tt), end=' ')
        #print('E %.4f' % (time.time()-T))
        #if modified:
        #    print('')
        #    print('end', id(self._group), id(self._pars))
        #print('group_sel:', self._group_sel)
        #print(mods_t['id'])
        #return True


    def _make_lines(self):
        self._lines_pref = self._lines_func.__name__+'_'+str(self._id)+'_'
        line = LMModel(self._lines_func, prefix=self._lines_pref,
                       series=self._series)
        d = self._defs
        self._pars = line.make_params()
        #print(d['z'])
        self._pars.add_many(
            #(self._lines_pref+'z', d['z'], d['z_vary'], 0, 10,
            (self._lines_pref+'z', d['z'], d['z_vary'], d['z']-d['z_min'], d['z']+d['z_max'],
             d['z_expr']),
            (self._lines_pref+'logN', d['logN'], d['logN_vary'], d['logN_min'],
             d['logN_max'], d['logN_expr']),
            (self._lines_pref+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._lines_pref+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']))
        self._lines = line


    def _make_lines_psf(self, spec, N_tot=False):
        if N_tot:
            lines_func = lines_voigt_N_tot
        else:
            lines_func = self._lines_func
        time_check = False
        if time_check:
            tt = time.time()
        self._lines_pref = self._lines_func.__name__+'_'+str(int(self._id))+'_'
        self._psf_pref = self._psf_func.__name__+'_'+str(int(self._id))+'_'
        if time_check:
            print('a %.4f' % (time.time()-tt))
            tt = time.time()
        line = LMModel(lines_func, prefix=self._lines_pref,
                       series=self._series)
        if time_check:
            print('b %.4f' % (time.time()-tt))
            tt = time.time()
        psf = LMModel(self._psf_func, prefix=self._psf_pref)
        if time_check:
            print('c %.4f' % (time.time()-tt))
            tt = time.time()
        line_psf = LMComposite(line, psf, convolve_simple)  # Time consuming
        if time_check:
            print('d %.4f' % (time.time()-tt))
            tt = time.time()

        d = self._defs

        if self._resol == None or np.isnan(self._resol):
            x = to_x(d['z'], trans_parse(self._series)[0])
            c = np.argmin(np.abs(spec.x.to(au.nm).value-x.to(au.nm).value))
            d['resol'] = spec.t['resol'][c]
            self._resol = d['resol']
        else:
            d['resol'] = self._resol

        self._pars = line_psf.make_params()
        self._pars.add_many(
            #(self._lines_pref+'z', d['z'], d['z_vary'], 0, 10,
            (self._lines_pref+'z', d['z'], d['z_vary'], d['z']-d['z_min'],
             d['z']+d['z_max'], d['z_expr']),
            (self._lines_pref+'logN', d['logN'], d['logN_vary'], d['logN_min'],
             d['logN_max'], d['logN_expr']),
            (self._lines_pref+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._lines_pref+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']),
            (self._psf_pref+'resol', d['resol'], d['resol_vary'],
             d['resol_min'], d['resol_max'], d['resol_expr']))

        if N_tot:
            self._pars.add_many(
                (self._lines_pref+'N_tot', 10**d['logN_min'], d['logN_vary'],
                 10**d['logN_min'], 10**d['logN_max'], ''),
                (self._lines_pref+'N_other', 10**d['logN_min'], True,
                 10**d['logN_min'], 10**d['logN_max'], ''))

        self._lines = line_psf

        if time_check:
            print('e %.4f' % (time.time()-tt))
            tt = time.time()

    def _make_regions(self, mod, spec, xs, thres, eval=False):
        #spec = mod._spec
        if 'fit_mask' not in spec.t.colnames:
            #logging.info("I'm adding column 'fit_mask' to spectrum.")
            spec.t['fit_mask'] = np.zeros(len(spec.x), dtype=bool)

        #self._pars.pretty_print()
        #print(xs)
        #print(mod)
        #print(self)
        #print(self._group)
        if eval:
            ys = mod.eval(x=xs, params=self._pars)
        else:
            ys = mod._ys
        #print('ys', len(ys))
        #print('xs', len(xs))
        c = []
        t = thres
        j = 0
        while len(c)==0 and j<1000:
            c = np.where(ys<1-t)[0]
            t = t*0.5
            j += 1

        #c = np.where(ys<1-thres)[0]
        #plt.plot(xs, ys)
        #plt.show()
        #if len(c)%2==0:
        #    c = c[:-1]
        xr = np.array(xs[c])
        yr = np.array(spec.y[c]/spec._t['cont'][c])
        #print(xr[np.argmin(yr)], yr[np.argmin(yr)])
        wr = np.array(spec._t['cont'][c]/spec.dy[c])
        spec.t['fit_mask'][c] = True
        #plt.plot(xr, mod.eval(x=xr, params=self._pars))
        #return xr, yr, wr, ys
        return xr, yr, wr



    def _make_N_tot(self, N_tot_specs=(None, None, None)):
        lines_pref, N_tot, N_other_expr = N_tot_specs
        if lines_pref is None: lines_pref = self._lines_pref
        pars = self._pars
        logN = pars[lines_pref+'logN']
        if N_tot is None:
            N_tot = np.sum([10**pars[p] for p in pars if 'logN' in p])
        N_other = N_tot-10**logN
        if N_other_expr is None:
            N_other_expr = ''
            for p in pars:
                if 'logN' in p and 'logN_' not in p and p != lines_pref+'logN':
                    N_other_expr += '+10**%s' % p
                    #N_other_expr += '%s' % p

        self._pars[lines_pref+'logN'].vary = False
        if lines_pref+'N_tot' not in self._pars:
            self._pars.add_many(
                (lines_pref+'N_tot', N_tot, True, N_other,
                 10**logN.max, ''))
        else:
            self._pars[lines_pref+'N_tot'].value = N_tot
            N_tot_min = N_other+0.5*10**self._pars[lines_pref+'logN']
            if N_tot_min < 1e10: N_tot_min = 1e10
            self._pars[lines_pref+'N_tot'].min = N_tot_min
        if lines_pref+'N_other' not in self._pars:
            self._pars.add_many(
                (lines_pref+'N_other', N_other, True, 10**logN.min,
                N_tot, N_other_expr))
        else:
            self._pars[lines_pref+'N_other'].value = N_other
            self._pars[lines_pref+'N_other'].expr = N_other_expr
            self._pars[lines_pref+'N_other'].max = N_tot\
                -0.5*10**self._pars[lines_pref+'logN']


    def _new_voigt(self, spec, series='Ly-a', z=2.0, logN=13, b=10, resol=None,
                   defs=None, N_tot=False, N_tot_specs=(None, None, None)):


        self._resol = resol
        self._series = series

        thres = defs.dict['group']['thres']

        time_check = False
        if time_check:
            tt = time.time()
        for l, v in zip(['z', 'logN', 'b', 'resol'], [z, logN, b, resol]):
            if l not in self._vars:
                self._vars[l] = v

        self._make_defs(defs.dict['voigt'])
        if time_check:
            print('a %.4f' % (time.time()-tt))
            tt = time.time()
        #print(self._defs)

        #N_tot = '47' == str(self._id)

        #self._make_lines()
        self._make_lines_psf(spec, N_tot)
        if time_check:
            print('b %.4f' % (time.time()-tt))
            tt = time.time()

        #self._make_group()
        self._make_groups(spec, thres=thres)
        if time_check:
            print('c %.4f' % (time.time()-tt))
            tt = time.time()

        self._make_comp2()
        if time_check:
            print('d %.4f' % (time.time()-tt))
            tt = time.time()

        ys = dc(self._ys)
        self._xr, self._yr, self._wr = \
            self._make_regions(self, spec, self._xs, thres=thres)
        if time_check:
            print('e %.4f' % (time.time()-tt))
            tt = time.time()
        self._xf, self._yf, self._wf = self._xr, self._yr, self._wr

        if N_tot:
            self._make_N_tot(N_tot_specs)

        #class_mute(self, Spectrum)
