from .functions import convolve, gaussian, voigt
from .vars import *
from astropy import units as au
from lmfit import CompositeModel as LMComposite
from lmfit import Model as LMModel
from lmfit import Parameters as LMParameters
#from lmfit.lineshapes import gaussian
import numpy as np

prefix = "Model:"

class Model(object):
    """ Class for models

    A Model is a combination of Lmfit Models for instrument PSF,
    continuum adjustment, and system profile."""

    def __init__(self,
                 spec):
        self._spec = spec

    def _create_comp(self, ion):
        lines = ModelLines(self._lines_func, ion=ion)
        psf = ModelPSF(self._psf_func)
        comp = LMComposite(lines, psf, convolve)
        comp._params = lines._params
        comp._params.update(psf._params)
        return comp
        #line = LMModel(getattr(f, self._line_func), prefix='line_', ion=ion)

    def single_voigt(self, series='Lya', z=2.0, N=1e14, b=20, btur=0,
                     resol=35000):
        """ @brief Create a voigt model of a single series
        @param series Series of transitions
        @param z Redshift
        @param N Column density
        @param b Doppler broadening
        @param btur Turbulence broadening
        @param resol Resolution
        @return 0
        """

        self._lines_func = voigt
        self._psf_func = gaussian
        self._comp = self._create_comp(series)

        x = np.array(self._spec._safe(self._spec.x).to(au.nm))
        chunk = np.where(np.logical_and(x > 413.5,
                                        x < 413.7))
        xc = np.array(x[chunk])# * self._spec.x.unit
        yc = np.array(self._spec.y[chunk]/self._spec._t['cont'][chunk])
        weights = np.array(self._spec._t['cont'][chunk]/self._spec.dy[chunk])

        self._comp._params.add_many(
            #('gaussian_c_min', int(0), False),
            #('gaussian_c_max', int(-1), False),
            ('gaussian_center', float(np.median(xc)), False),
            ('gaussian_resol', float(resol), False),
            ('voigt_z', float(z)),
            ('voigt_N', float(N)),
            ('voigt_b', float(b)),
            ('voigt_btur', float(btur)))
        self._comp._params.pretty_print()
        guess = self._comp.eval(self._comp._params, x=x)
        print(len(guess), np.sum(self._spec._where_safe))
        self._spec.t['model'] = self._spec.y
        self._spec.t['model'][self._spec._where_safe] = guess * self._spec._t['cont'][self._spec._where_safe]
        print("guess done")
        xc = np.array(self._spec.x[chunk].to(au.nm))# * self._spec.x.unit
        yc = np.array(self._spec.y[chunk]/self._spec._t['cont'][chunk])
        weights = np.array(self._spec._t['cont'][chunk]/self._spec.dy[chunk])
        fit = self._comp.fit(yc, self._comp._params, x=xc)
        print("fit done")
        self._spec.t['model'][chunk] = fit.best_fit * self._spec._t['cont'][chunk]
        print(fit.fit_report())
        return 0

class ModelLines(LMModel):
    """ Class for line models

    A ModelLines is a model for spectral line."""

    def __init__(self,
                 func,
                 **kwargs):

        if func.__name__ != 'voigt':
            print(prefix, "Only 'voigt' function supported for lines.")
            return None
        prefix = func.__name__+'_'
        super(ModelLines, self).__init__(func, prefix=prefix, **kwargs)
        self._params = self.make_params()

class ModelPSF(LMModel):
    """ Class for psf models

    A ModelPSF is a model for the instrumental PSF."""

    def __init__(self,
                 func,
                 **kwargs):

        if func.__name__ != 'gaussian':
            print(prefix, "Only 'gaussian' function supported for lines.")
            return None
        prefix = func.__name__+'_'
        super(ModelPSF, self).__init__(func, prefix=prefix, **kwargs)
        self._params = self.make_params()
