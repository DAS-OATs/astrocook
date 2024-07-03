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


    def clip_flux(self, ran='all', smooth_len=200, kappa_hi=6, kappa_lo=3,
                  iter=100, fudge=1, delta_x=5000):
        """ @brief Clip flux
        @details Discriminate absorbed spectrum bins by applying a kappa-sigma
        clipping within a running window.

        The recipes computes the mean of `y` in a running window across the
        spectrum and saves it in column `y_rm`. It then computes the deviation
        $$\Delta y =$$ `y` $$-$$ `y_rm` for all spectral bins and finds the
        outliers in this distribution.

        Since the recipe is meant to operate on quasar spectra, it assumes
        that the absorption features are narrow and widespread, while the
        emission feature are wide enough to be considered a part of the
        continuum. For this reason, the recipe is biased towards detecting
        outliers with a *negative* deviation (lower outliers, i.e. absorbed
        bins) instead of a *positive* deviation (upper outliers, most likely
        residuals of sky line or cosmic ray hit subtraction).

        In practice, a bin is regarded as an outlier if it has either
        $$\Delta y < -$$`kappa_hi` $$\\times$$ `dy` or
        $$\Delta y > $$`kappa_lo` $$\\times$$ `dy`.

        The clipping is iterated several times, always removing the outliers
        and re-computing `y_rm`, until no more outliers are found or a maximum
        number of iterations `iter` is reached.

        A column `y_abs` is created and set equal to `1` for all the lower
        outliers (absorbed bins), `0` elsewhere. Similarly, a column `y_em` is
        created for the upper outliers (“emitted” bins) and a column `y_cont`
        is created for the remaining bins (continuum bins). The latter is
        multiplied for a `fudge` factor and smoothed in the velocity space with
        a gaussian filter with standard deviation `std`. The procedure is
        applied in the range `xmin`-`xmax`, except for the smoothing which is
        applied to the whole spectral range. The result is saved in column
        `cont` of the spectrum as the current estimate of the emission
        continuum. A set of `delta_x`-spaced nodes is superimposed to the final
        continuum
        @param ran Wavelength range (nm)
        @param smooth_len Smoothing length (km/s)
        @param kappa_hi Number of standard deviations for clipping above
        @param kappa_lo Number of standard deviations for clipping below
        @param fudge Fudge factor to scale the continuum
        @param delta_x Spacing of nodes (km/s)
        @return 0
        """

        try:
            xmin, xmax = parse_range(ran)
            smooth_len = float(smooth_len)
            kappa_hi = float(kappa_hi)
            kappa_lo = float(kappa_lo)
            fudge = float(fudge)
            delta_x = float(delta_x)
        except:
            logging.error(msg_param_fail)
            return 0

        maxiter = 1000

        spec = self.sess.spec
        dv = spec._dv()

        # Prepare columns
        if 'mask_abs' not in spec._t.colnames:
            spec._t['mask_abs'] = at.Column(np.zeros(len(spec._t)), dtype=int)
        if 'cont' not in spec._t.colnames:
            spec._t['cont'] = at.Column(np.array(None, ndmin=1), dtype=float)


        # Extract range
        xsel = np.logical_and(spec._t['x']>xmin, spec._t['x']<xmax)
        xsels = np.where(xsel==1)[0][0]
        spec_t = spec._t[xsel]

        # Compute running median
        hwindow = int(smooth_len/np.min(dv[xsel]))//8
        y_rm = running_mean(spec_t['y'], h=hwindow)

        # Clip flux
        sum_sel = len(spec_t)
        sel = np.zeros(len(spec_t), dtype=bool)
        for i in range(maxiter):

            sel[~sel] = spec_t['y'][~sel]-y_rm<-kappa_lo*spec_t['dy'][~sel]
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

        spec._t['mask_abs']
        cont_temp = np.interp(spec._t['x'][xsel], x_rm, y_rm)*spec_t['y'].unit

        # Apply fudge
        cont_temp *= fudge

        # Merge with old continuum
        if 'cont' in spec._t.colnames:
            cont_old = spec._t['cont']
        else:
            cont_old = [np.nan]*len(spec._t)
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

        return 0
