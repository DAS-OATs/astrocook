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
                 pars=None):
        self._series = series
        self._lines_func = lines_func
        self._psf_func = psf_func
        #self._cont_func = cont_func
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        self._pars = pars
        self._set_defaults()
        self._create_comp()

    def _create_comp(self):
        #cont = ModelCont(self._cont_func)
        lines = ModelLines(self, **self._lines_d)#._lines_func, series=self._series)
        psf = ModelPSF(self, **self._psf_d)
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


    def single_voigt(self, series='Ly_a', z=2.0, N=1e14, b=20, btur=0,
                     resol=35000, cont_ampl=0.1):
        """ @brief Create a voigt model of a single series
        @param series Series of transitions
        @param z Redshift
        @param N Column density
        @param b Doppler broadening
        @param btur Turbulence broadening
        @param resol Resolution
        @return 0
        """

        z = float(z)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)
        cont_ampl = float(cont_ampl)

        spec = self._sess.spec

        #self._cont_func = cont_gauss
        self._lines_func = voigt
        self._psf_func = gaussian
        self._comp = self._create_comp(series)

        trans = series_d[series]

        x = np.array(spec._safe(spec.x).to(au.nm))
        s = spec._where_safe
        c = np.array([], dtype=int)
        #"""
        for t in series_d[series]:
            c = np.append(c, np.where(np.logical_and(x > (0.995+z)*xem_d[t].value,
                                                     x < (1.005+z)*xem_d[t].value))[0])
        """
        c = np.where(np.logical_and(x > (0.995+z)*xem_d[trans[0]].value,
                                    x < (1.005+z)*xem_d[trans[-1]].value))
        """
        xc = np.array(x[c])# * spec.x.unit
        if 'noabs' in spec._t.colnames:
            yc = np.array(spec._t['noabs'][c]/spec._t['cont'][c])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        w = np.array(spec._t['cont'][c]/spec.dy[c])

        self._comp._pars.add_many(
            #('gaussian_center', float(np.median(xc)), False),
            #('cont_gauss_amplitude', cont_ampl, True),
            #('cont_gauss_sigma', 100, True),
            ('lines_voigt_z', z, True),
            ('lines_voigt_N', N, True),
            ('lines_voigt_b', b, True),
            ('lines_voigt_btur', btur, True),
            ('psf_gaussian_resol', resol, False),
            ('psf_gaussian_z', z, False))

        guess = self._comp.eval(self._comp._pars, x=x)

        if 'model' in spec._t.colnames:
            print(prefix, "I'm updating column 'model'.")
        else:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*spec.y.unit
            spec._t['model'][s] = spec._t['cont'][s]



        fit = self._comp.fit(yc, self._comp._pars, x=xc)
        spec._t['model'][c] = fit.best_fit * spec._t['model'][c]#spec._t['cont'][c]
        if 'noabs' in spec._t.colnames:
            print(prefix, "I'm updating column 'noabs'.")
        else:
            print(prefix, "I'm adding column 'noabs'.")
            spec._t['noabs'] = spec.y
            spec._t['noabs'][c] = spec._t['cont'][c]+spec.y[c]-spec._t['model'][c]
        return 0

class ModelCont(LMModel):
    """ Class for continuum models

    A ModelCont is a model for the local continuum."""

    def __init__(self,
                 func,
                 **kwargs):

        if func.__name__ != 'gaussian':
            print(prefix, "Only 'gaussian' function supported for lines.")
            return None
        prefix = 'cont_'+func.__name__+'_'
        super(ModelCont, self).__init__(func, prefix=prefix, **kwargs)
        self._pars = self.make_pars()

class ModelLines(LMModel):
    """ Class for line models

    A ModelLines is a model for spectral line."""

    def __init__(self, call,
                 **kwargs):
        func = call._lines_func
        func_name = func.__name__
        if func_name != 'lines_voigt':
            print(prefix, "Only 'lines_voigt' function supported for lines.")
            return None
        prefix = func_name+'_'
        super(ModelLines, self).__init__(func, prefix=prefix)#, **kwargs)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_lines_voigt(self, call):
        d = call._lines_d
        self._pars = self.make_params()
        self._pars.add_many(
            ('lines_voigt_z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            ('lines_voigt_N', d['N'], d['N_vary'], d['N_min'], d['N_max'],
             d['N_expr']),
            ('lines_voigt_b', d['b'], d['b_vary'], d['b_min'], d['b_max'],
             d['b_expr']),
            ('lines_voigt_btur', d['btur'], d['btur_vary'], d['btur_min'],
             d['btur_max'], d['btur_expr']))

class ModelPSF(LMModel):
    """ Class for psf models

    A ModelPSF is a model for the instrumental PSF."""

    def __init__(self, call,
                 **kwargs):
        func = call._psf_func
        func_name = func.__name__
        if func_name != 'psf_gauss':
            print(prefix, "Only 'psf_gauss' function supported for PSF.")
            return None
        prefix = func_name+'_'
        super(ModelPSF, self).__init__(func, prefix=prefix)#, **kwargs)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_psf_gauss(self, call):
        d = call._psf_d
        self._pars = self.make_params()
        self._pars.add_many(
            ('psf_gauss_z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            ('psf_gauss_resol', d['resol'], d['resol_vary'], d['resol_min'],
             d['resol_max'], d['resol_expr']))
