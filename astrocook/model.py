from .functions import convolve, psf_gauss, lines_voigt
from .vars import *
from astropy import units as au
from astropy.io import ascii
from lmfit import CompositeModel as LMComposite
from lmfit import Model as LMModel
from lmfit import Parameters as LMParameters
#from lmfit.lineshapes import gaussian as gauss
import numpy as np

prefix = "Model:"

class Model(object):
    """ Class for models

    A Model is a combination of Lmfit Models for instrument PSF,
    continuum adjustment, and system profile."""

    def __init__(self,
                 series='Ly_a',
                 lines_func=lines_voigt,
                 psf_func=psf_gauss,
                 #cont_func=cont_gauss,
                 z=None,
                 zmin=None,
                 zmax=None,
                 pars=None,
                 count=0):
        self._series = series
        self._lines_func = lines_func
        self._psf_func = psf_func
        self._cont_func = cont_func
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        self._pars = pars
        self._count = count
        self._set_defaults()
        self._create_comp()

    def _create_comp(self):
        cont = ModelCont(self)
        lines = ModelLines(self, self._count, self._series, **self._lines_d)#._lines_func, series=self._series)
        psf = ModelPSF(self, self._count, **self._psf_d)
        #cont_lines = cont * lines
        #comp = LMComposite(cont_lines, psf, convolve)
        self._comp = LMComposite(lines, psf, convolve)
        self._comp._pars = lines._pars
        #comp._pars.update(cont._pars)
        self._comp._pars.update(psf._pars)

    def _set_defaults(self):

        # Line defaults
        self._lines_d = lines_voigt_d
        self._lines_d['z'] = self._z
        self._lines_d['z_min'] = self._zmin
        self._lines_d['z_max'] = self._zmax
        self._psf_d = psf_gauss_d
        self._psf_d['z'] = self._z
        self._psf_d['zmin'] = self._zmin
        self._psf_d['zmax'] = self._zmax
        if self._pars is not None:
            for p in self._pars:
                if p in self._lines_d:
                    self._lines_d[p] = self._pars[p]
                if p in self._psf_d:
                    self._psf_d[p] = self._pars[p]


class ModelCont(LMModel):
    """ Class for continuum models

    A ModelCont is a model for the local continuum."""

    def __init__(self, call, count,
                 **kwargs):
        func = call._cont_func
        func_name = func.__name__
        if func_name != 'cont_gauss':
            print(prefix, "Only 'cont_gauss' function supported for lines.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelCont, self).__init__(func, prefix=self._prefix)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_cont_gauss(self, call):
        d = call._lines_d
        self._pars = self.make_params()
        self._pars.add_many(
            (self._prefix+'ampl', d['ampl']))

class ModelLines(LMModel):
    """ Class for line models

    A ModelLines is a model for spectral line."""

    def __init__(self, call, count, series,
                 **kwargs):
        func = call._lines_func
        func_name = func.__name__
        if func_name != 'lines_voigt':
            print(prefix, "Only 'lines_voigt' function supported for lines.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelLines, self).__init__(func, prefix=self._prefix,
                                         series=series)#, **kwargs)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_lines_voigt(self, call):
        d = call._lines_d
        self._pars = self.make_params()
        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'N', d['N'], d['N_vary'], d['N_min'], d['N_max'],
             d['N_expr']),
            (self._prefix+'b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            (self._prefix+'btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']))

class ModelPSF(LMModel):
    """ Class for psf models

    A ModelPSF is a model for the instrumental PSF."""

    def __init__(self, call, count,
                 **kwargs):
        func = call._psf_func
        func_name = func.__name__
        if func_name != 'psf_gauss':
            print(prefix, "Only 'psf_gauss' function supported for PSF.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelPSF, self).__init__(func, prefix=self._prefix)#, **kwargs)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_psf_gauss(self, call):
        d = call._psf_d
        self._pars = self.make_params()
        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'resol', d['resol'], d['resol_vary'], d['resol_min'],
             d['resol_max'], d['resol_expr']))
