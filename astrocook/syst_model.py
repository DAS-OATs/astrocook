from .functions import adj_gauss, lines_voigt, convolve, psf_gauss
from .vars import *
from lmfit import CompositeModel as LMComposite
from lmfit import Model as LMModel
from lmfit import Parameters as LMParameters
import numpy as np

prefix = "System model:"


class SystModel(LMComposite):
    """ Class for composite models

    A composite model is a combination of Lmfit Models for instrument PSF,
    continuum adjustment, and system profile."""
    def __init__(self, spec, series, vars, regs=[],
                 lines_func=lines_voigt,
                 psf_func=psf_gauss,
                 cont_func=None,
                 **kwargs):

        self._spec = spec
        self._series = series
        self._vars = vars
        self._lines_func = lines_func
        self._psf_func = psf_func
        self._make_defs()
        self._make_lines()
        self._make_regs()
        self._make_psf()
        self._make_comp()
        self._pars.pretty_print()


        #self._make_pars(

    def _make_comp(self):
        self._group = np.prod(self._lines)
        super(SystModel, self).__init__(self._group, self._psf, convolve)

    def _make_defs(self):
        self._defs = pars_std_d
        for v in self._vars:
            if v in self._defs:
                self._defs[v] = self._vars[v]

    def _make_lines(self):
        self._lines = []
        self._lines_pref = self._lines_func.__name__+'_0_'
        line = LMModel(self._lines_func, prefix=self._lines_pref,
                       series=self._series)
        d = self._defs
        self._pars = line.make_params()
        self._pars.add_many(
            (self._lines_pref+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._lines_pref+'logN', d['logN'], d['logN_vary'], d['logN_min'],
             d['logN_max'], d['logN_expr']),
            (self._lines_pref+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._lines_pref+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']))
        self._lines.append(line)

    def _make_psf(self):
        for i, r in enumerate(self._xr):
            self._psf_pref = self._psf_func.__name__+'_'+str(i)+'_'
            if i == 0:
                self._psf = LMModel(self._psf_func, prefix=self._psf_pref,
                                    reg=r)
            else:
                self._psf += LMModel(self._psf_func, prefix=self._psf_pref,
                                     reg=r)
            d = self._defs
            self._pars.update(self._psf.make_params())
            self._pars.add_many(
                (self._psf_pref+'resol', d['resol'], d['resol_vary'],
                 d['resol_min'], d['resol_max'], d['resol_expr']))

    def _make_regs(self):
        self._xr = []
        self._yr = []

        """
        spec = self._spec
        xs = np.array(spec._safe(spec.x).to(au.nm))
        if 'deabs' in spec._t.colnames:
            ys = np.array(spec._t['deabs']/spec._t['cont'])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        wc = np.array(spec._t['cont'][c]/spec.dy[c])
        for i, r in enumerate([self._xr]):
        """

        self._xr.append(np.array(self._spec._safe(self._spec.x).to(au.nm)))
        self._yr.append(np.array(self._spec._safe(self._spec.x).to(au.nm)))

    def fit(self):
        xs = np.ravel(self._xr)
