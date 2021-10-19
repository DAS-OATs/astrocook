from .functions import *
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

thres = 1e-3

class SystModel(LMComposite):

    def __init__(self, spec, systs, series=[], vars=None, constr=None, z0=None,
                 lines_func=lines_voigt,
                 psf_func=psf_gauss,
                 cont_func=None):
        self._spec = spec
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


    def _fit(self, fit_kws={}):
        vary = np.any([self._pars[p].vary for p in self._pars])
        #print(vary)
        if vary:
            time_start = datetime.datetime.now()
            #for p in self._pars:
            #    if '_z' in p:
            #        print(self._pars[p])
            #self._pars.pretty_print()
            if 'max_nfev' in fit_kws:
                max_nfev = fit_kws['max_nfev']
                del fit_kws['max_nfev']
            else:
                max_nfev = None
            fit_kws['ftol'] = 1e-3
            #print('out', len(self._xf), self._xf)
            fit = super(SystModel, self).fit(self._yf, self._pars, x=self._xf,
                                             weights=self._wf,
                                             max_nfev=max_nfev,
                                             fit_kws=fit_kws,
                                             nan_policy='omit',
                                             #fit_kws={'method':'lm'},
                                             method='least_squares')
                                             #method='emcee')
            #print(fit.result.success)
            #print(fit.result.message)
            time_end = datetime.datetime.now()
            self._pars = fit.params
            #for p in self._pars:
            #    if '_z' in p:
            #        print(self._pars[p])
            #self._pars.pretty_print()
            self._ys = self.eval(x=self._xs, params=self._pars)
            self._chi2r = fit.redchi
            self._aic = fit.aic
            self._bic = fit.bic
            return 0
        else:
            return 1


    def _make_comp(self):
        super(SystModel, self).__init__(self._group, self._psf, convolve_simple)
        #self._pars.pretty_print()

    def _zero(self, x):
        return 0*x

    def _make_comp2(self):
        super(SystModel, self).__init__(self._group, LMModel(self._zero), operator.add)
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

    def _make_group(self, thres=thres):
        """ @brief Group lines that must be fitted together into a single model.
        """

        spec = self._spec

        mods_t = self._mods_t
        self._xs = np.array(spec._safe(spec.x).to(au.nm))
        #ys = self._lines.eval(x=self._xs, params=self._pars)

        ys = self._lines.eval(x=self._xs, params=self._pars)
        #self._pars.pretty_print()

        self._group = self._lines
        self._group_list = []
        #plt.plot(self._xs, ys)

        for i, s in enumerate(mods_t):
            mod = s['mod']
            #ys_s = mod.eval(x=self._xs, params=mod._pars)
            ys_s = mod._ys
            ymax = np.maximum(ys, ys_s)
            #print(np.amin(ymax))
            y_cond = np.amin(ymax)<1-thres or np.amin(ymax)==np.amin(ys)
            pars_cond = False
            #for p,v in self._pars.items():
            for p,v in self._constr.items():
                for mod_p,mod_v in mod._pars.items():
                    pars_cond = pars_cond or v==mod_p
            #print(s['id'],pars_cond)
            if y_cond or pars_cond:
                """
                print('mod')
                mod._pars.pretty_print()
                print('self')
                self._pars.pretty_print()
                """
                self._group *= mod._group
                mod = s['mod']
                for p,v in mod._pars.items():
                    if v.expr != None:
                        self._constr[p] = v.expr
                        v.expr = ''

                """
                print('mod cleaned')
                mod._pars.pretty_print()
                """
                self._pars.update(mod._pars)
                """
                print('combine')
                self._pars.pretty_print()
                """
                #self._pars.pretty_print()
                #print(self._constr)
                if pars_cond or self._constr != {}:
                    for p,v in self._constr.items():
                        #print(s['id'], pars_cond, self._constr != {}, self._constr, p, v)
                        #print(self._pars[p])
                        #print(self._pars[v])
                        self._pars[p].expr = v
                        if v != '':
                            #print(self._pars[v])
                            try:
                                self._pars[p].min = self._pars[v].min
                                self._pars[p].max = self._pars[v].max
                                self._pars[p].value = self._pars[v].value
                            except:
                                self._pars[p].expr = ''
                        #print(self._pars[p])
                """
                print('constrained')
                self._pars.pretty_print()
                """
                self._group_list.append(i)
                #self._pars.pretty_print()

                ##
                mod._ys = self._group.eval(x=self._xs, params=self._pars)
                ##
        #print(self._group_list)
        #print('final')
        #self._pars.pretty_print()
        if len(self._group_list) > 1:
            ids = [i for il in mods_t['id'][self._group_list[1:]] for i in il]
            mods_t.remove_rows(self._group_list[1:])
            for i in ids:
                mods_t['id'][self._group_list[0]].append(i)
        if self._group_list == []:
            self._group_sel = -1
        else:
            self._group_sel = self._group_list[0]
        self._ys = self._group.eval(x=self._xs, params=self._pars)
        #plt.plot(self._xs, self._ys, linestyle='--')
        #plt.show()


    def _make_group2(self, thres=thres):
        """ @brief Group lines that must be fitted together into a single model.
        """
        spec = self._spec

        mods_t = self._mods_t
        self._xs = np.array(spec._safe(spec.x).to(au.nm))

        #print(self.__dict__)
        ys = self._lines.eval(x=self._xs, params=self._pars)
        ysmin = np.amin(ys)
        self._group = self._lines
        self._group_list = []

        #print(mods_t['id'])
        modified = False
        for i, s in enumerate(mods_t):
            ttt = time.time()
            mod = s['mod']
            ys_s = mod._ys

            #print('0', time.time()-ttt)
            ttt = time.time()
            cond, pars_cond = False, False
            for p,v in self._constr.items():
                for mod_p,mod_v in mod._pars.items():
                    cond = cond or v==mod_p
            if cond: pars_cond = True
            #print('1', time.time()-ttt)
            #ttt = time.time()
            if not cond:
            #print(s['id'])
                #ymax = np.maximum(ys, ys_s)
                yminmax = np.amin(np.maximum(ys,ys_s))
                #print('2', time.time()-ttt)
                #ttt = time.time()
                #cond = np.amin(ymax)<1-thres or np.amin(ymax)==np.amin(ys)
                cond = yminmax<1-thres or yminmax==ysmin
                #print('3', time.time()-ttt)
                #ttt = time.time()

            #print('1', time.time()-ttt)
            #ttt = time.time()
            if cond: #y_cond or pars_cond:
                #print('mod')
                #mod._pars.pretty_print()
                #print('self')
                #self._pars.pretty_print()
                #try:
                self._group *= mod._group
                #except:
                #    self._group_sel = -1
                #    return
                mod = s['mod']
                for p,v in mod._pars.items():
                    if v.expr != None:
                        self._constr[p] = v.expr
                        v.expr = ''
                self._pars.update(mod._pars)
                #print('2', time.time()-ttt)
                #ttt = time.time()
                if pars_cond or self._constr != {}:
                    for p,v in self._constr.items():
                        self._pars[p].expr = v
                        if v != '':
                            try:
                                self._pars[p].min = self._pars[v].min
                                self._pars[p].max = self._pars[v].max
                                self._pars[p].value = self._pars[v].value
                            except:
                                self._pars[p].expr = ''
                self._group_list.append(i)
                #print('3', time.time()-ttt)
                #ttt = time.time()
                #mod._ys = self._group.eval(x=self._xs, params=self._pars)
                mod._ys = ys_s*ys
                modified = True
                #print('4', time.time()-ttt)
                #ttt = time.time()
                #print('')
                #print(i, id(self._group), id(self._pars))

        if len(self._group_list) > 1:
            ids = [i for il in mods_t['id'][self._group_list[1:]] for i in il]
            mods_t.remove_rows(self._group_list[1:])
            for i in ids:
                mods_t['id'][self._group_list[0]].append(i)
        if self._group_list == []:
            self._group_sel = -1
        else:
            self._group_sel = self._group_list[0]
        if modified and False:
            self._ys = mod._ys
            print('hey')
        else:
            self._ys = self._group.eval(x=self._xs, params=self._pars)
        #if modified:
        #    print('')
        #    print('end', id(self._group), id(self._pars))


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


    def _make_lines_psf(self):
        self._lines_pref = self._lines_func.__name__+'_'+str(self._id)+'_'
        self._psf_pref = self._psf_func.__name__+'_'+str(self._id)+'_'
        line = LMModel(self._lines_func, prefix=self._lines_pref,
                       series=self._series)
        psf = LMModel(self._psf_func, prefix=self._psf_pref, spec=self._spec)
        line_psf = LMComposite(line, psf, convolve_simple)

        d = self._defs

        if self._resol == None or np.isnan(self._resol):
            #c = np.where(self._spec.x.to(au.nm).value==self._xs[len(self._xs)//2])
            #d['resol'] = self._spec.t['resol'][c][0]
            x = to_x(d['z'], trans_parse(self._series)[0])
            c = np.argmin(np.abs(self._spec.x.to(au.nm).value-x.to(au.nm).value))
            d['resol'] = self._spec.t['resol'][c]
        else:
            d['resol'] = self._resol

        self._pars = line_psf.make_params()
        #print(d['z'])
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

        self._lines = line_psf


    def _make_psf(self):
        d = self._defs
        """
        for i, r in enumerate(self._xr):
            if self._resol == None:
                c = np.where(self._spec.x.to(au.nm).value==r[len(r)//2])
                d['resol'] = self._spec.t['resol'][c][0]
            self._psf_pref = self._psf_func.__name__+'_'+str(i)+'_'
            psf = LMModel(self._psf_func, prefix=self._psf_pref, reg=r)
            if i == 0:
                self._psf = psf
            else:
                self._psf += psf
            self._pars.update(psf.make_params())
            self._pars.add_many(
                (self._psf_pref+'resol', d['resol'], d['resol_vary'],
                 d['resol_min'], d['resol_max'], d['resol_expr']))
        """

        if self._resol == None:
            #c = np.where(self._spec.x.to(au.nm).value==self._xs[len(self._xs)//2])
            #d['resol'] = self._spec.t['resol'][c][0]
            x = to_x(d['z'], trans_parse(self._series)[0])
            c = np.argmin(np.abs(self._spec.x.to(au.nm).value-x.to(au.nm).value))
            d['resol'] = self._spec.t['resol'][c]
        else:
            d['resol'] = self._resol

        self._psf_pref = self._psf_func.__name__+'_0_'
        psf = LMModel(self._psf_func, prefix=self._psf_pref, spec=self._spec)
        self._psf = psf
        self._pars.update(psf.make_params())
        self._pars.add_many(
            (self._psf_pref+'resol', d['resol'], d['resol_vary'],
             d['resol_min'], d['resol_max'], d['resol_expr']))
        #"""

    def _make_regions(self, mod, xs, thres=thres):
        spec = self._spec
        if 'fit_mask' not in spec.t.colnames:
            #logging.info("I'm adding column 'fit_mask' to spectrum.")
            spec.t['fit_mask'] = np.zeros(len(spec.x), dtype=bool)


        #self._pars.pretty_print()
        #print(xs)
        #print(mod)
        #print(self)
        #print(self._group)
        #ys = mod.eval(x=xs, params=self._pars)
        ys = self._ys
        c = []
        t = thres
        while len(c)==0:
            c = np.where(ys<1-t)[0]
            t = t*0.5

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



    def _make_regs(self, thres=thres):
        spec = self._spec

        #ys = self._group.eval(x=self._xs, params=self._pars)
        #c = np.where(ys<1-thres)[0]
        c = []
        t = thres
        while len(c)==0:
            c = np.where(self._ys<1-t)[0]
            t = t*0.5

        self._xr = np.split(self._xs[c], np.where(np.ediff1d(c)>1.5)[0]+1)

        self._xf = np.concatenate([np.array(x) for  x in self._xr])
        self._yf = np.array(spec.y[c]/spec._t['cont'][c])
        #self._wf = np.array(spec._t['cont'][c]/spec.dy[c])
        self._wf = np.array(np.power(spec._t['cont'][c]/spec.dy[c], 2))
        try:
            self._xm = np.concatenate([np.arange(x[0], x[-1], 1e-5) \
                                      for x in self._xr])
        except:
            self._xm = np.array([])

    def _new_voigt(self, series='Ly-a', z=2.0, logN=13, b=10, resol=None, defs=None):
        #if resol == None:
        #    self._resol = self._spec.t['resol'][len(self._spec.t)//2]
        #else:
        self._resol = resol
        self._series = series
        tt = time.time()
        for l, v in zip(['z', 'logN', 'b', 'resol'], [z, logN, b, resol]):
            if l not in self._vars:
                self._vars[l] = v

        self._make_defs(defs)
        #print('a', time.time()-tt)
        #tt = time.time()

        #self._make_lines()
        self._make_lines_psf()
        #print('b', time.time()-tt)
        #tt = time.time()

        #self._make_group()
        self._make_group2()
        #print('c', time.time()-tt)
        #tt = time.time()

        #self._make_psf()
        #self._make_comp()
        self._make_comp2()
        #print('d', time.time()-tt)
        #tt = time.time()
        ys = dc(self._ys)
        #self._xr, self._yr, self._wr, self._ys = self._make_regions(self, self._xs)
        self._xr, self._yr, self._wr = self._make_regions(self, self._xs)
        if np.sum(self._ys-ys)>0: print('attenzione:', np.sum(self._ys-ys))
        #print('e', time.time()-tt)
        #tt = time.time()

        #print('r', self._yr)
        self._xf, self._yf, self._wf = self._xr, self._yr, self._wr
        #self._pars.pretty_print()
        #self._ys = self.eval(x=self._xs, params=self._pars)
