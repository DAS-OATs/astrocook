from . import Spec1D, Line
from .utils import *
from astropy import units as u
from astropy.table import Column, Table
from copy import deepcopy as dc
import numpy as np
from scipy.stats import sigmaclip
from statsmodels.nonparametric.smoothers_lowess import lowess

class Cont(Spec1D, Line):
    def __init__(self, spec=None, line=None,
                 x=None, y=None, dy=None,
                 xunit=xunit_def, yunit=yunit_def,
                 meta=None, dtype=float):  
        """ Constructor for the Cont class """ 

        self._spec = spec
        self._line = line
        if (x is not None):
            if (xunit == None):
                x = x * xunit_def
        if (y is not None):
            if (yunit == None):
                y = y * yunit_def
        if (dy is not None):      
            if (yunit == None):
                dy = dy * yunit_def
        if (x is not None and y is not None):
            self._t = self.create_t(x, y, dy)

        self._use_good = False
        
    def create_t(self, x, y, dy=None,
                 xunit=xunit_def, yunit=yunit_def, dtype=float):
        """ Create a spectrum of the continuum """

        x = np.array(x, ndmin=1)
        y = np.array(y, ndmin=1)
        dy = np.array(dy, ndmin=1)
        t = Table()
        t['X'] = Column(x, dtype=dtype, unit=xunit)
        t['Y'] = Column(y, dtype=dtype, unit=yunit)
        t['DY'] = Column(dy, dtype=dtype, unit=yunit)

        return t

    def smooth_maxima(self, smooth=3.0, flux_corr=1.0, kappa_low=2.0,
                    kappa_high=2.0):
        """ Determine the emission continuum by removing absorption lines """
        
        range_x = np.max(self._spec.x) - np.min(self._spec.x)
        x = self._line._maxs['X']
        y = self._line._maxs['Y'] * flux_corr

        clip_x = x
        clip_y = y
        stop = False
        i = 0

        while (stop == False):
            frac = smooth*u.nm/range_x 
            le = lowess(clip_y, clip_x, frac=frac, it=0, delta=0.0,
                        is_sorted=True)
            cont_y = np.interp(clip_x, le[:, 0], le[:, 1])
            norm_y = clip_y / cont_y
            clip_y = sigmaclip(norm_y, low=kappa_low, high=kappa_high)[0]
            clip_x = clip_x[np.isin(norm_y, clip_y)]
            cont_y = cont_y[np.isin(norm_y, clip_y)]
            clip_y = clip_y * cont_y
            stop = ((len(clip_y) == len(norm_y)) \
                    or (len(clip_y) < 100))

        x = self._spec.x
        y = np.interp(x, le[:, 0], le[:, 1]) * self._spec.y.unit

        self._t = self.create_t(x, y)
