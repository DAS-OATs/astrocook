from .functions import adj_gauss, lines_voigt, convolve, psf_gauss
from .vars import *
from astropy import table as at
from copy import deepcopy as dc
import datetime
from lmfit import CompositeModel as LMComposite
from lmfit import Model as LMModel
from lmfit import Parameters as LMParameters
from matplotlib import pyplot as plt
import numpy as np

thres = 5e-3

class SystModel(LMComposite):

    def __init__(self, spec, systs, series=[], vars=[], z0=None,
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
        self._vars = vars
        self._z0 = z0
        self._lines_func = lines_func
        self._psf_func = psf_func


    def _fit(self, fit_kws={}):
        time_start = datetime.datetime.now()
        fit = super(SystModel, self).fit(self._yf, self._pars, x=self._xf,
                                         weights=self._wf, fit_kws=fit_kws,
                                         method='least_squares')
                                         #method='emcee')
        time_end = datetime.datetime.now()
        #print(fit.nfev, time_end-time_start)
        self._pars = fit.params
        self._chi2r = fit.redchi
        self._aic = fit.aic
        self._bic = fit.bic

    def _make_comp(self):
        super(SystModel, self).__init__(self._group, self._psf, convolve)
        #self._pars.pretty_print()

    def _make_defs(self):
        self._defs = pars_std_d
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

        try:
            ok
    #    except:
            d = self._defs
            psf = LMModel(self._psf_func, prefix='temp_', reg=self._xs)
            temp = dc(self)
            temp._pars.update(psf.make_params())
            temp._pars.add_many(('temp_resol', d['resol'], d['resol_vary'],
                                d['resol_min'], d['resol_max'], d['resol_expr']))
            #temp._pars.pretty_print()
            comp = LMComposite(temp._lines, psf, convolve)
            ys = comp.eval(x=self._xs, params=temp._pars)
        except:
            ys = self._lines.eval(x=self._xs, params=self._pars)
        #self._pars.pretty_print()

        self._group = self._lines
        self._group_list = []
        #plt.plot(self._xs, ys)
        for i, s in enumerate(mods_t):
            mod = s['mod']
            #mod._pars.pretty_print()
            ys_s = mod.eval(x=self._xs, params=mod._pars)
            #ys_s = mod._ys
            #print(s['z0'], s['id'])
            #print(np.amin(np.maximum(ys, ys_s)) < 1-thres)
            ymax = np.maximum(ys, ys_s)
            #plt.plot(self._xs, ys, linestyle='--')
            #plt.plot(self._xs, ys_s, linestyle=':')
            #plt.plot(self._xs, ymax)
            #plt.scatter(self._xs[np.argmin(ymax)], np.min(ymax))
            #plt.show()
            if np.amin(ymax)<1-thres or np.amin(ymax)==np.amin(ys):
                self._group *= mod._group
                self._pars.update(mod._pars)
                self._group_list.append(i)
                #self._pars.pretty_print()

        if len(self._group_list) > 1:
            #print(self._group_list)
            #print(mods_t['id'][self._group_list[1:]])
            #ids = list(np.ravel([np.array(i)
            #                     for i in mods_t['id'][self._group_list[1:]]]))
            ids = [i for il in mods_t['id'][self._group_list[1:]] for i in il]
            #print(ids)
            mods_t.remove_rows(self._group_list[1:])
            for i in ids:
                #mods_t['id'][0].append(i)
                mods_t['id'][self._group_list[0]].append(i)
        if self._group_list == []:
            self._group_sel = -1
        else:
            self._group_sel = self._group_list[0]
        self._ys = self._group.eval(x=self._xs, params=self._pars)
        #plt.plot(self._xs, self._ys)
        #plt.show()

    def _make_lines(self):
        self._lines_pref = self._lines_func.__name__+'_'+str(self._id)+'_'
        line = LMModel(self._lines_func, prefix=self._lines_pref,
                       series=self._series)
        d = self._defs
        self._pars = line.make_params()
        self._pars.add_many(
            #(self._lines_pref+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
            # d['z_expr']),
            (self._lines_pref+'z', d['z'], d['z_vary'], d['z']-1e-4, d['z']+1e-4,
             d['z_expr']),
            (self._lines_pref+'logN', d['logN'], d['logN_vary'], d['logN_min'],
             d['logN_max'], d['logN_expr']),
            (self._lines_pref+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._lines_pref+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']))
        self._lines = line


    def _make_psf(self):
        d = self._defs
        for i, r in enumerate(self._xr):
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
        self._wf = np.array(spec._t['cont'][c]/spec.dy[c])
        try:
            self._xm = np.concatenate([np.arange(x[0], x[-1], 1e-5) \
                                      for x in self._xr])
        except:
            self._xm = np.array([])

    def _new_voigt(self, series='Ly_a', z=2.0, logN=13, b=10, resol=70000):
        self._series = series
        self._vars = {'z': z, 'logN': logN, 'b': b, 'resol': resol}
        self._make_defs()
        self._make_lines()
        self._make_group()
        self._make_regs()
        self._make_psf()
        self._make_comp()
        #self._ys = self.eval(x=self._xs, params=self._pars)
