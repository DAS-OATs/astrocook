from .functions import running_mean
from .message import *
from .utils import parse_range

import astropy.constants as ac
import astropy.table as at
import astropy.units as au
import logging

class CookbookContinuum(object):

    def __init__(self):
        super(CookbookContinuum, self).__init__()


    def clip_flux(self, ran='all', smooth_len=400, kappa=2, fudge='auto',
                  knots_dist=2000, mode='update'):
        """ @brief Clip flux
        @details Estimate the continuum by clipping absorbers.

        @url continuum_cb.html#clip-flux

        @param ran Wavelength range (nm)
        @param smooth_len Smoothing length (km/s)
        @param kappa Number of sigma to reject absorber
        @param fudge Fudge factor to scale the continuum
        @param knots_dist Distance between knots (km/s)
        @param mode Update or replace
        @return 0
        """

        try:
            xmin, xmax = parse_range(ran)
            smooth_len = float(smooth_len)
            fudge = None if fudge=='auto' else float(fudge)
            kappa = float(kappa)
            knots_dist = float(knots_dist)
            mode = str(mode)
        except:
            logging.error(msg_param_fail)
            return 0

        if mode not in ['update', 'replace']:
            logging.warning("I cannot understand the mode. I will use `update`.")

        maxiter = 1000

        spec = self.sess.spec
        dv = spec._dv()

        # Extract range
        xsel = np.logical_and(spec._t['x']>xmin, spec._t['x']<xmax)
        xsels = np.where(xsel==1)[0][0]
        spec_t = spec._t[xsel]

        # Prepare columns
        if 'mask_abs' not in spec._t.colnames:
            logging.info("I'm adding column `mask_abs`.")
        if 'cont' not in spec._t.colnames:
            logging.info("I'm adding column `cont`.")
        if mode=='replace' or 'mask_abs' not in spec._t.colnames:
            spec._t['mask_abs'] = at.Column(np.zeros(len(spec._t)), dtype=int)
        else:
            spec._t['mask_abs'][xsel] = 0
        if mode=='replace' or 'cont' not in spec._t.colnames:
            spec._t['cont'] = at.Column(np.array(None, ndmin=1), dtype=float)


        # Compute running median
        hwindow = int(smooth_len/np.min(dv[xsel]))//8
        y_rm = running_mean(spec_t['y'], h=hwindow)

        # Clip flux
        sum_sel = len(spec_t)
        sel = np.zeros(len(spec_t), dtype=bool)
        for i in range(maxiter):

            sel[~sel] = spec_t['y'][~sel]-y_rm<-kappa*spec_t['dy'][~sel]
            selw = np.where(sel==1)[0]+xsels
            spec._t['mask_abs'][selw] = np.ones(len(selw))
            x_rm = spec_t['x'][~sel]
            y_rm = running_mean(spec_t['y'][~sel], h=hwindow)

            if i == maxiter-1 and sum_sel != np.sum(sel):
                logging.warning("Clipping not converged after {} iterations! "\
                                .format(i+1))
            if sum_sel != np.sum(sel):
                sum_sel = np.sum(sel)
            else:
                logging.info("Clipping converged after %i iterations." % (i+1))
                break

        # Compute fudge if `auto` and apply it
        if fudge is None:
            diff = running_mean(spec_t['y'][~sel]-y_rm, hwindow)
            fudge = 1+diff/y_rm
        y_rm *= fudge

        cont_temp = np.interp(spec._t['x'][xsel], x_rm, y_rm)*spec_t['y'].unit

        # Update old continuum
        if mode=='update' and 'cont' in spec._t.colnames:
            cont_old = spec._t['cont']
        else:
            cont_old = np.array([np.nan]*len(spec._t))
        cont_old[xsel] = cont_temp
        spec._t['cont_old'] = cont_old
        spec._gauss_convolve(std=smooth_len, input_col='cont_old',
                             output_col='cont')
        spec._t.remove_column('cont_old')

        # Propagate to lines
        lines = self.sess.lines
        if lines is not None:
            lines._t['cont'] = np.interp(lines._t['x'], spec._t['x'],
                                         spec._t['cont'])

        # Extract nodes
        self.nodes_extract(delta_x=knots_dist, mode='cont')

        return 0
