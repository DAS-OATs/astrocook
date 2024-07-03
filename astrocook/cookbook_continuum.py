class CookbookContinuum(object):

    def __init__(self):
        super(CookbookContinuum, self).__init__()


    def clip_flux(self, xmin, xmax, hwindow=100, kappa_hi=6, kappa_lo=3,
                  iter=100, fudge=1, std=500, delta_x=5000):
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
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @param hwindow Half-window size in pixels for running mean
        @param kappa_hi Number of standard deviations for clipping above
        @param kappa_lo Number of standard deviations for clipping below
        @param iter Number of iterations
        @param fudge Fudge factor to scale the continuum
        @param std Standard deviation for gaussian convolution (km/s)
        @param delta_x Spacing of nodes (km/s)
        @return 0
        """

        try:
            xmin = 0*au.nm if xmin in ["", "None", None] else float(xmin) * au.nm
            xmax = 1e4*au.nm if xmax in ["", "None", None] else float(xmax) * au.nm
            hwindow = int(hwindow)
            kappa_hi = float(kappa_hi)
            kappa_lo = float(kappa_lo)
            iter = int(iter)
            fudge = float(fudge)
            std = float(std)
            delta_x = float(delta_x)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        lines = self.sess.lines

        xsel = np.logical_and(spec._t['x']>xmin, spec._t['x']<xmax)
        spec_t = spec._t[xsel]

        y_rm = running_mean(spec_t['y'], h=hwindow)

        if 'y_rm' not in spec._t.colnames:
            logging.info("I'm adding column 'y_rm'.")
        else:
            logging.warning("I'm updating column 'y_rm'.")
        spec_t['y_rm'] = at.Column(y_rm, dtype=float)

        if 'y_em' not in spec_t.colnames:
            logging.info("I'm adding column 'y_em'.")
        spec_t['y_em'] = at.Column(np.ones(len(spec_t)), dtype=int)
        if 'y_abs' not in spec_t.colnames:
            logging.info("I'm adding column 'y_abs'.")
        spec_t['y_abs'] = at.Column(np.zeros(len(spec_t)), dtype=int)
        if 'y_cont' not in spec_t.colnames:
            logging.info("I'm adding column 'y_cont'.")
        spec_t['y_cont'] = spec_t['y']
        if 'cont' not in spec_t.colnames:
            logging.info("I'm adding column 'cont'.")
        spec_t['cont'] = at.Column(np.array(None, ndmin=1), dtype=float)
        #spec_t['cont'].unit = spec_t['y'].unit

        sum_sel = len(spec_t)
        for i in range(iter):
            sel = spec_t['y']-spec_t['y_rm']<-kappa_lo*spec_t['dy']
            sel = np.logical_or(sel, spec_t['y']-spec_t['y_rm']>kappa_hi*spec_t['dy'])
            spec_t['y_em'][sel] = 0
            spec_t['y_abs'][sel] = 1
            x_rm = spec_t['x'][~sel]
            if 2*hwindow+1>len(spec_t['y'][~sel]):
                logging.warning("Too many outliers to compute a running "
                                "median at iteration %i! Try increasing sigma."
                                % i)
                y_rm = np.ones(len(spec_t['y'][~sel])) \
                       * np.median(spec_t['y'][~sel])
            else:
                y_rm = running_mean(spec_t['y'][~sel], h=hwindow)

            if i == iter-1 and sum_sel != np.sum(sel):
                logging.warning("Clipping not converged after %i iterations! "
                                "Try iterating more." % (i+1))
            if sum_sel != np.sum(sel):
                sum_sel = np.sum(sel)
            else:
                logging.info("Clipping converged after %i iterations." % (i+1))
                break

            spec_t['y_rm'] = np.interp(spec_t['x'], x_rm, y_rm)
        x_cont = spec_t['x'][np.where(spec_t['y_em'])]
        y_cont = spec_t['y'][np.where(spec_t['y_em'])]
        spec_t['y_cont'] = np.interp(spec_t['x'], x_rm, y_rm)*spec_t['y'].unit
        if lines is not None:
            lines._t['y_cont'] = np.interp(lines._t['x'], x_rm, y_rm)*lines._t['y'].unit

        spec_t['y_cont'] *= fudge
        if 'cont' in spec._t.colnames:
            spec._t['y_cont'] = spec._t['cont']
        else:
            spec._t['y_cont'] = np.nan
        spec._t['y_cont'][xsel] = spec_t['y_cont']

        spec._gauss_convolve(std=std, input_col='y_cont', output_col='cont')
        self.nodes_extract(delta_x=delta_x, mode='cont')

        return 0
