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

    def spec_from_struct(self, x='0,spec,x', y='0,spec,y', dy='0,spec,y'):
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
        @param dy Expression for the error on y
        @return Session with synthetic spectrum
        """

        for i, s in enumerate(self.sess._gui._sess_list):
            if s.spec is not None:
                for c in sorted(s.spec._t.colnames, key=len, reverse=True):
                    xold, yold = dc(x), dc(y)
                    x = x.replace('%i,spec,%s' % (i, c),
                                  str(list(np.array(s.spec._t[c]))))
                    y = y.replace('%i,spec,%s' % (i, c),
                                  str(list(np.array(s.spec._t[c]))))
                    dy = dy.replace('%i,spec,%s' % (i, c),
                                    str(list(np.array(s.spec._t[c]))))
                    if x != xold:
                        xunit = s.spec._t[c].unit
                    if y != yold:
                        yunit = s.spec._t[c].unit

        x = expr_eval(ast.parse(x, mode='eval').body)
        y = expr_eval(ast.parse(y, mode='eval').body)
        dy = expr_eval(ast.parse(dy, mode='eval').body)

        xmin, xmax = create_xmin_xmax(x)
        #dy = y/snr

        # Add gaussian noise
        rng = np.random.default_rng()
        norm = rng.standard_normal(size=y.size)
        y = y+dy*norm

        spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit)
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_synth',
                      spec=spec)

        return new

    def spec_from_systs(self, x='0,spec,x', dy='0,spec,dy', sess='0',
                        resol=None):
        """@brief Create a synthetic spectrum from systems
        @details Create a synthetic spectrum from a list of systems taken from
        an existing session
        @param x Expression for wavelength-like array
        @param dy Expression for the error on y
        @param sess Number of the session with the systems
        @param resol Resolution
        @return Session with synthetic spectrum
        """

        if resol is not None:
            resol = float(resol)

        for i, s in enumerate(self.sess._gui._sess_list):
            if s.spec is not None:
                for c in sorted(s.spec._t.colnames, key=len, reverse=True):
                    xold, dyold = dc(x), dc(dy)
                    x = x.replace('%i,spec,%s' % (i, c),
                                  str(list(np.array(s.spec._t[c]))))
                    dy = dy.replace('%i,spec,%s' % (i, c),
                                    str(list(np.array(s.spec._t[c]))))
                    if x != xold:
                        xunit = s.spec._t[c].unit
                    if dy != dyold:
                        yunit = s.spec._t[c].unit

        x = expr_eval(ast.parse(x, mode='eval').body)
        dy = expr_eval(ast.parse(dy, mode='eval').body)
        xmin, xmax = create_xmin_xmax(x)

        parse = self.sess._gui._panel_sess._struct_parse(sess+',systs')
        if parse is None: return 0
        _, systs, _ = parse

        y = np.ones(len(x))
        for i, r in enumerate(systs._mods_t):
            mod = r['mod']
            if resol is not None:
                mod._pars['psf_gauss_%i_resol' % mod._id].value = resol
            #mod._pars.pretty_print()
            y = mod.eval(x=x, params=mod._pars) * y

        # Add gaussian noise
        rng = np.random.default_rng()
        norm = rng.standard_normal(size=y.size)
        y = y+dy*norm
        #print(len(x), len(y), len(dy), len(xmin), len(xmax))

        spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit)
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_synth',
                      spec=spec)

        new._systs = systs

        return new
