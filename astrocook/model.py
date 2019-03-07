from .functions import adj_gauss, convolve, lines_voigt, psf_gauss
from .vars import *
from astropy import units as au
from astropy.io import ascii
from lmfit import CompositeModel as LMComposite
from lmfit import Model as LMModel
from lmfit import Parameters as LMParameters
#from lmfit.lineshapes import gaussian as gauss
import numpy as np

prefix = "Model:"

class Model(LMModel):
    """ Class for models

    A Model is a combination of Lmfit Models for instrument PSF,
    continuum adjustment, and system profile."""

    def __init__(self, func, count, series, z, zmin, zmax):
        self._func = func
        self._count = count
        self._series = series
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        prefix = func.__name__+'_'+str(count)+'_'
        super(Model, self).__init__(func, prefix=prefix, series=self._series)

    def _make_pars(self, func, **kwargs):

        # Line defaults
        d = pars_d[func.__name__+'_d']
        d['z'] = self._z
        d['z_min'] = self._zmin
        d['z_max'] = self._zmax
        if kwargs is not None:
            for p in kwargs:
                if p in d:
                    d[p] = kwargs[p]
        self._pars = self.make_params()
        self._pars.pretty_print()
        getattr(self, '_make_pars_'+func.__name__)(d, **kwargs)

    def _make_pars_lines_voigt(self, d, **kwargs):

        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'N', d['N'], d['N_vary'], d['N_min'], d['N_max'],
             d['N_expr']),
            (self._prefix+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._prefix+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']))

    def _make_pars_psf_gauss(self, d, **kwargs):

        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'resol', d['resol'], d['resol_vary'], d['resol_min'],
             d['resol_max'], d['resol_expr']))

class ModelLines(Model):
    """ Class for line models

    A ModelLines is a model for spectral line."""

    def __init__(self, func, count, series, z, zmin, zmax, **kwargs):
        if func.__name__ != 'lines_voigt':
            print(prefix, "Only 'lines_voigt' function supported for lines.")
            return None
        super(ModelLines, self).__init__(func, count, series, z, zmin, zmax)
        self._make_pars(func, **kwargs)



class ModelPSF(Model):
    """ Class for PSF models

    A ModelLines is a model for the instrument PSF."""

    def __init__(self, func, count, series, z, zmin, zmax, **kwargs):
        if func.__name__ != 'psf_gauss':
            print(prefix, "Only 'psf_gauss' function supported for lines.")
            return None
        super(ModelPSF, self).__init__(func, count, series, z, zmin, zmax)
        self._make_pars(func, **kwargs)
