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

class Model(object):
    """ Class for models

    A Model is a combination of Lmfit Models for instrument PSF,
    continuum adjustment, and system profile."""

    def __init__(self,
                 series='Ly_a',
                 lines_func=lines_voigt,
                 psf_func=psf_gauss,
                 adj_func=adj_gauss,
                 z=None,
                 zmin=None,
                 zmax=None,
                 pars=None,
                 count=0):
        self._series = series
        self._lines_func = lines_func
        self._psf_func = psf_func
        self._adj_func = adj_func
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        self._pars = pars
        self._count = count
        self._set_defaults()
        self._create_comp()

    def _create_comp(self):
        #adj = ModelAdj(self, self._count, **self._adj_d)
        lines = ModelLines(self, self._count, self._series, **self._lines_d)#._lines_func, series=self._series)
        psf = ModelPSF(self, self._count, **self._psf_d)
        #adj_lines = adj * lines
        adj_lines = lines
        self._comp = LMComposite(adj_lines, psf, convolve)
        #self._comp._pars = adj._pars
        self._comp._pars = lines._pars
        #self._comp._pars.update(adj._pars)
        self._comp._pars.update(psf._pars)

    def _set_defaults(self):

        # Line defaults
        """
        self._adj_d = adj_gauss_d
        self._adj_d['z'] = self._z
        self._adj_d['z_min'] = self._zmin
        self._adj_d['z_max'] = self._zmax
        """
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
                """
                if p in self._adj_d:
                    self._adj_d[p] = self._pars[p]
                """
                if p in self._lines_d:
                    self._lines_d[p] = self._pars[p]
                if p in self._psf_d:
                    self._psf_d[p] = self._pars[p]


class ModelAdj(LMModel):
    """ Class for continuum models

    A ModelAdj is a model for the local continuum."""

    def __init__(self, call, count,
                 **kwargs):
        func = call._adj_func
        func_name = func.__name__
        if func_name != 'adj_gauss':
            print(prefix, "Only 'adj_gauss' function supported for lines.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelAdj, self).__init__(func, prefix=self._prefix)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_adj_gauss(self, call):
        d = call._adj_d
        self._pars = self.make_params()
        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'ampl', d['ampl'], d['ampl_vary'], d['ampl_min'],
             d['ampl_max'], d['ampl_expr']),
            (self._prefix+'sigma', d['sigma'], d['sigma_vary'], d['sigma_min'],
             d['sigma_max'], d['sigma_expr']))

class ModelLines(LMModel):
    """ Class for line models

    A ModelLines is a model for spectral line."""

    #def __init__(self, call, count, series,
    #             **kwargs):
    #    func = call._lines_func
    def __init__(self, func, series, z, zmin, zmax, count, **kwargs):
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        self._count = count
        func_name = func.__name__
        if func_name != 'lines_voigt':
            print(prefix, "Only 'lines_voigt' function supported for lines.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelLines, self).__init__(func, prefix=self._prefix,
                                         series=series)#, **kwargs)
        #getattr(self, '_pars_'+func_name)(call)
        getattr(self, '_pars_'+func_name)(**kwargs)

    #def _pars_lines_voigt(self, call):
    #    d = call._lines_d
    def _pars_lines_voigt(self, **kwargs):

        # Line defaults
        d = lines_voigt_d
        d['z'] = self._z
        d['z_min'] = self._zmin
        d['z_max'] = self._zmax
        if kwargs is not None:
            for p in kwargs:
                if p in d:
                    d[p] = kwargs[p]
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
