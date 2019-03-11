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
    def __init__(self, count, series, z, zmin, zmax,
                 lines_func=lines_voigt,
                 psf_func=psf_gauss,
                 cont_func=None,
                 **kwargs):

        self._count = count
        self._series = series
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        self._lines_func = lines_func
        self._psf_func = psf_func
        self._psf_pref = psf_func.__name__+'_'+str(self._count)+'_'
        self._lines_pref = lines_func.__name__+'_'+str(self._count)+'_'
        self._new()
        self._comp()

    def _comp(self):
        self._group = np.prod(self._lines)
        super(SystModel, self).__init__(self._group, self._psf, convolve)

    def _new(self):
        lines_pref = self._lines_func.__name__+'_'+str(self._count)+'_'
        psf_pref = self._psf_func.__name__+'_'+str(self._count)+'_'
        self._lines = [LMModel(self._lines_func, prefix=self._lines_pref,
                               series=self._series)]

        self._psf = LMModel(self._psf_func, prefix=self._psf_pref,
                            series=self._series)

    def add(self):
        self._count += 1


class SystModelStd(SystModel):

    def __init__(self, count, series, z, zmin, zmax, **kwargs):

        super(SystModelStd, self).__init__(count, series, z, zmin, zmax,
                                           **kwargs)
        self._make_pars(**kwargs)
        #self._pars.pretty_print()
        #print(self._pars)


    def _make_pars(self, **kwargs):

        # Line defaults
        print(kwargs)
        d = pars_std_d
        d['z'] = self._z
        d['z_min'] = self._zmin
        d['z_max'] = self._zmax
        if kwargs is not None:
            for p in kwargs:
                if p in d:
                    d[p] = kwargs[p]
        self._pars = self.make_params()
        self._pars.add_many(
            (self._lines_pref+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._lines_pref+'N', d['N'], d['N_vary'], d['N_min'], d['N_max'],
             d['N_expr']),
            (self._lines_pref+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._lines_pref+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']),
            (self._psf_pref+'z', d['z'], False, d['z_min'], d['z_max'],
             d['z_expr']),
            (self._psf_pref+'resol', d['resol'], d['resol_vary'], d['resol_min'],
             d['resol_max'], d['resol_expr']))
