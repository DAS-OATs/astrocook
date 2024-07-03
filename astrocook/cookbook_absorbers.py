from .vars import resol_def
from .line_list import LineList

import logging
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import argrelmin, argrelmax, find_peaks
from scipy.special import erf, erfc, erfinv
import sys
import time

prefix = "[INFO] cookbook_absorbers:"

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        super(CookbookAbsorbers, self).__init__()
        self._refit_n = 3
        self._chi2rav_thres = 1e-2
        self._chi2r_thres = np.inf
        self._dlogN_thres = np.inf
        self._max_nfev = max_nfev_def
        self._sel_fit = False

    def _lines_cands_find(self, series, z_start, z_end, dz):
        return self.sess.lines._cands_find(series, z_start, z_end, dz)


    def find_lines(self, kind='abs', prominence=None, append=True):
        """ @brief Find lines

        @details Find absorption or emission lines, based on their prominence.

        @param kind Kind of lines (`abs` or `em`)
        @param prominence Prominence of lines (as in `scipy.signal.find_peaks`)
        @param append Append lines to existing line list
        @return 0
        """

        try:
            kind = str(kind)
            prominence = None if prominence in [None, 'None'] \
                else float(prominence)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if kind not in ['abs', 'em']:
            logging.error("`kind` should be `abs` or `em`. Aborting.")
            return 0

        spec = self.sess.spec
        fact = -1 if kind=='abs' else 1

        ynorm = fact*(spec._t['y'])
        if prominence is None: prominence = 5*(spec._t['dy'])

        peaks, properties = find_peaks(ynorm, prominence=prominence)
        lines = LineList(row=spec._t[peaks], source='y', kind=kind,
                         xunit=spec._xunit, yunit=spec._yunit, meta=spec._meta)
        lines.append_replace(append, self.sess)

        return 0
