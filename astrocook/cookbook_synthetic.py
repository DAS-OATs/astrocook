from .functions import create_xmin_xmax, expr_eval
from .message import *
from .spectrum import Spectrum
import ast
from copy import deepcopy as dc
import logging
import numpy as np


class CookbookSynthetic(object):
    """ Cookbook of utilities to create and handle synthetic spectra
    """

    def __init__(self):
        super(CookbookSynthetic, self).__init__()

    def spec_from_struct(self, x='0,spec,x', y='0,spec,y', snr=100, ron=0.01):
        """@brief Create a synthetic spectrum
        @details Create a synthetic spectrum from existing structures (a
        wavelenght-like array and a flux-like array). The structure expressions
        must be parsable by AST, with columns described by a string with the
        session number, the structure tag (spec, lines, systs), and the column
        name separated by a comma (e.g. 0,spec,x, meaning "column x of spectrum
        from session 0"). A gaussian noise is added to the spectrum to match a
        given signal-to-noise ratio. A new session is created with the synthetic
        spectrum.
        @param x Expression for wavelength-like array
        @param y Expression for flux-like array
        @param snr Signal-to-noise ratio (single value or expression)
        @param ron Read-out noise
        @param ks Kernel size
        @return Session with synthetic spectrum
        """

        try:
            snr = float(snr)
            snr_expr = False
        except:
            snr_expr = True

        for i, s in enumerate(self.sess._gui._sess_list):
            if s.spec is not None:
                for c in sorted(s.spec._t.colnames, key=len, reverse=True):
                    xold, yold = dc(x), dc(y)
                    x = x.replace('%i,spec,%s' % (i, c),
                                  str(list(np.array(s.spec._t[c]))))
                    y = y.replace('%i,spec,%s' % (i, c),
                                  str(list(np.array(s.spec._t[c]))))
                    if snr_expr:
                        snr = snr.replace('%i,spec,%s' % (i, c),
                                          str(list(np.array(s.spec._t[c]))))
                    if x != xold:
                        xunit = s.spec._t[c].unit
                    if y != xold:
                        yunit = s.spec._t[c].unit

        x = expr_eval(ast.parse(x, mode='eval').body)
        y = expr_eval(ast.parse(y, mode='eval').body)
        if snr_expr: snr = expr_eval(ast.parse(snr, mode='eval').body)

        xmin, xmax = create_xmin_xmax(x)
        dy = y/snr

        # Add gaussian noise
        rng = np.random.default_rng()
        norm = rng.standard_normal(size=y.size)
        y = y+np.sqrt(dy**2+ron**2)*norm

        spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit)
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_synth',
                      spec=spec)

        return new
