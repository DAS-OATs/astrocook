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
                 sess):
        self._sess = sess

    def _create_comp(self, series):
        lines = ModelLines(self._lines_func, series=series)
        psf = ModelPSF(self._psf_func)
        comp = LMComposite(lines, psf, convolve)
        comp._params = lines._params
        comp._params.update(psf._params)
        return comp
        #line = LMModel(getattr(f, self._line_func), prefix='line_', ion=ion)

    def single_voigt(self, series='Ly_a', z=2.0, N=1e14, b=20, btur=0,
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

        z = float(z)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)

        spec = self._sess.spec

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

        self._comp._params.add_many(
            #('gaussian_center', float(np.median(xc)), False),
            ('gaussian_resol', resol, False),
            ('gaussian_z', z, False),
            ('voigt_z', z, True),
            ('voigt_N', N, True),
            ('voigt_b', b, True),
            ('voigt_btur', btur, True))

        guess = self._comp.eval(self._comp._params, x=x)

        if 'model' in spec._t.colnames:
            print(prefix, "I'm updating column 'model'.")
        else:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*spec.y.unit
            spec._t['model'][s] = spec._t['cont'][s]



        fit = self._comp.fit(yc, self._comp._params, x=xc)
        spec._t['model'][c] = fit.best_fit * spec._t['model'][c]#spec._t['cont'][c]
        if 'noabs' in spec._t.colnames:
            print(prefix, "I'm updating column 'noabs'.")
        else:
            print(prefix, "I'm adding column 'noabs'.")
            spec._t['noabs'] = spec.y
            spec._t['noabs'][c] = spec._t['cont'][c]+spec.y[c]-spec._t['model'][c]
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
