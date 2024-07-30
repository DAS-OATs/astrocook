from .functions import create_xmin_xmax, expr_eval, resol_check
from .message import *
from .spectrum import Spectrum
from .vars import *
import ast
from copy import deepcopy as dc
import logging
import numpy as np


class CookbookSynthetic(object):
    """ Cookbook of utilities to create and handle synthetic spectra
    """

    def __init__(self):
        super(CookbookSynthetic, self).__init__()

    def spec_from_struct(self, x='0,spec,x', y='0,spec,y', dy='0,spec,dy'):
        """@brief Synthetic spectrum from structures
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
        """@brief Synthetic spectrum from systems
        @details Create a synthetic spectrum from a list of systems taken from
        an existing session.
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

        #parse = self.sess._gui._panel_sess._struct_parse(sess+',systs')
        parse = self._struct_parse(sess+',systs')
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


    def spec_from_systs_random(self, n, series='Ly-a',
                               z_min=0, z_max=6, z_seed=None,
                               logN_mean=12.5,
                               logN_std=0.7, logN_seed=None,
                               b_min=pars_std_d['b_min'],
                               b_max=pars_std_d['b_max'], b_seed=None,
                               resol=resol_def):
        """ @brief Synthetic spectrum from random systems
        @details Create a synthetic spectrum from a list of systems with random
        redshifts, column density, and Doppler broadening.
        @param n Number of systems
        @param series Series of transitions
        @param z_min Minimum redshift
        @param z_max Maximum redshift
        @param z_seed Seed for random sampling in [z_min, z_max]
        @param logN_mean Mean (logarithmic) column density
        @param logN_std Standard deviation of (logarithmic) column density
        @param logN_seed Seed for random sampling in [logN_min, logN_max]
        @param b_min Minimum Doppler broadening
        @param b_max Maximum Doppler broadening
        @param b_seed Seed for random sampling in [b_min, b_max]
        @param resol Resolution
        @return Session with synthetic spectrum
        """

        try:
            n = int(n)
            z_min = float(z_min)
            z_max = float(z_max)
            z_seed = None if z_seed in [None, 'None'] else int(z_seed)
            logN_mean = float(logN_mean)
            logN_std = float(logN_std)
            logN_seed = None if logN_seed in [None, 'None'] else int(logN_seed)
            b_min = float(b_min)
            b_max = float(b_max)
            b_seed = None if b_seed in [None, 'None'] else int(b_seed)
            resol = None if resol in [None, 'None'] else float(resol)
            """
            if snr in [None, 'None']:
                snr = None
            elif snr in [np.inf, 'inf']:
                snr = np.inf
            else:
                snr = float(snr)
            append = str(append) == 'True'
            """
        except:
            logging.error(msg_param_fail)
            return 0

        #self._chi2r_thres = float(chi2r_thres)
        #self._dlogN_thres = float(dlogN_thres)
        self._refit_n = 0
        #self._chi2rav_thres = float(chi2rav_thres)
        self._max_nfev = 0 #float(max_nfev)


        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0

        for s in series.split(';'):
            if z_max != z_min:
                rng = np.random.default_rng(z_seed)
                z_list = rng.random(n)*(z_max-z_min)+z_min
            else:
                z_list = np.full(n, z_min)
            """
            if logN_max != logN_min:
                rng = np.random.default_rng(logN_seed)
                logN_list = rng.random(n)*(logN_max-logN_min)+logN_min
                logN_list = rng.normal(logN_mean, logN_std, n)
            else:
                logN_list = np.full(n, logN_min)
            """
            rng = np.random.default_rng(logN_seed)
            logN_list = rng.normal(logN_mean, logN_std, n)
            if b_max != b_min:
                rng = np.random.default_rng(b_seed)
                b_list = rng.random(n)*(b_max-b_min)+b_min
            else:
                b_list = np.full(n, b_min)

            s_list = [s]*n
            resol_list = [resol]*n


            self._systs_prepare(False)
            self._systs_add(s_list, z_list, logN_list, b_list,
                            resol_list=resol_list)
            self._systs_cycle()
            self._spec_update()

        spec = self.sess.spec
        x, xmin, xmax = spec._t['x'], spec._t['xmin'], spec._t['xmax']
        """
        if snr is None:
            rng = np.random.default_rng()
            norm = rng.standard_normal(size=x.size)
            dy = spec._t['dy']
            y = (spec._t['model']+dy*norm)/spec._t['cont']*spec._t['y']
        else:
            y = spec._t['model']
            dy = y/snr
        if 'cont' in spec._t.colnames:
            y, dy = y/spec._t['cont'], dy/spec._t['cont']
        xunit = spec._t['x'].unit
        yunit = spec._t['y'].unit
        spec_new = Spectrum(x, xmin, xmax, y, dy, xunit, yunit)
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_synth',
                      spec=spec_new)

        new._systs = self.sess.systs

        return new
        """

        spec._t['y'] *= spec._t['model']/spec._t['cont']
        return 0
