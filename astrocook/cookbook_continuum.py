from .functions import *
from .functions import _voigt_par_convert, _fadd
from .message import *
from .vars import *
from astropy import table as at
from astropy import units as au
#from scipy.interpolate import UnivariateSpline as uspline
from scipy.optimize import root_scalar
from matplotlib import pyplot as plt

class CookbookContinuum(object):
    """ Cookbook of utilities for continuum fitting
    """

    def __init__(self):
        super(CookbookContinuum, self).__init__()


    def _voigt_fwhm(self, l):
        x = np.arange(l['x']-1.0, l['x']+1.0, 1e-4)#*au.nm
        mod = 1-lines_voigt(x, l['z'], l['logN'], l['b'], 0.0, l['series'])
        hm = 0.5*(np.max(mod))
        xw = x[np.where(mod>hm)]
        fwhm = xw[-1]-xw[0]
        return fwhm


    def node_add(self, x, y):
        self.sess.spec._node_add(self.sess.nodes, x, y)


    def node_remove(self, x):
        self.sess.spec._node_remove(self.sess.nodes, x)


    def node_interp(self):
        self.sess.spec._node_interp(self.sess.nodes, self.sess.lines)


### Basic

    def nodes_clean(self, kappa=5.0):
        """ @brief Clean nodes
        @details Clean the list of nodes from outliers.
        @param kappa Number of standard deviation away from the window average
        @return 0
        """
        try:
            kappa = float(kappa)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.nodes = self.sess.spec._nodes_clean(self.sess.nodes, kappa)
        return 0

    def nodes_extract(self, delta_x=500, xunit=au.km/au.s):
        """ @brief Extract nodes
        @details Extract nodes from a spectrum. Nodes are averages of x and y in
        slices, computed after masking lines.
        @param delta_x Size of slices
        @param xunit Unit of wavelength or velocity
        @return 0
        """
        try:
            xunit = au.Unit(xunit)
            delta_x = float(delta_x)*xunit
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.nodes = self.sess.spec._nodes_extract(delta_x, xunit)
        #print(len(self.sess.nodes.t))

        return 0


    def nodes_interp(self, smooth=0):
        """ @brief Interpolate nodes
        @details Interpolate nodes with a univariate spline to estimate the
        emission level.
        @param smooth Smoothing of the spline
        @return 0
        """

        try:
            smooth = float(smooth)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._nodes_interp(self.sess.lines, self.sess.nodes)
        return 0


    def peaks_find(self, col='conv', kind='min', kappa=5.0, append=True):
        """ @brief Find peaks
        @details Find the peaks in a spectrum column. Peaks are the extrema
        (minima or maxima) that are more prominent than a given number of
        standard deviations. They are saved as a list of lines.
        @param col Column where to look for peaks
        @param kind Kind of extrema ('min' or 'max')
        @param kappa Number of standard deviations
        @param append Append peaks to existing line list
        @return 0
        """

        try:
            kappa = float(kappa)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        self._peaks_found = False
        spec = self.sess.spec
        if col not in spec.t.colnames:
            logging.error("The spectrum has not a column named '%s'. Please "\
                          "pick another one." % col)
            return 0

        peaks = spec._peaks_find(col, kind, kappa)
        if len(peaks.t) > 0:
            self._peaks_found = True
            source = [col]*len(peaks.t)
            from .line_list import LineList
            lines = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.t[col],
                             peaks.dy, source, spec._xunit, spec._yunit,
                             meta=spec._meta)
            if append and self.sess.lines is not None \
                and len(self.sess.lines.t) > 0:
                self.sess.lines._append(lines)
                self.sess.lines._clean()
            else:
                self.sess.lines = lines

            self.sess.lines_kind = 'peaks'
            spec._lines_mask(self.sess.lines, source=col)
        return 0


### Advanced


    def flux_clip(self, hwindow=100, kappa=5, iter=100, std=500):
        """ @brief Clip flux in spectrum
        @details Select flux pixels to find outliers
        @param hwindow Half-window size in pixels for running mean
        @param kappa Number of standar deviations for clipping
        @param iter Number of iterations
        @param std Standard deviation for gaussian convolution (km/s)
        @return 0
        """

        try:
            hwindow = int(hwindow)
            kappa = float(kappa)
            iter = int(iter)
            std = float(std)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        lines = self.sess.lines

        y_rm = running_mean(spec._t['y'], h=hwindow)
        #plt.plot(spec._t['x'], spec._t['y'])
        #plt.plot(spec._t['x'], y_rm)
        #plt.show()

        if 'y_rm' not in spec._t.colnames:
            logging.info("I'm adding column 'y_rm'.")
        else:
            logging.warning("I'm updating column 'y_rm'.")
        spec._t['y_rm'] = at.Column(y_rm, dtype=float)

        if 'y_em' not in spec._t.colnames:
            logging.info("I'm adding column 'y_em'.")
        spec._t['y_em'] = at.Column(np.ones(len(spec._t)), dtype=int)
        if 'y_abs' not in spec._t.colnames:
            logging.info("I'm adding column 'y_abs'.")
        spec._t['y_abs'] = at.Column(np.zeros(len(spec._t)), dtype=int)
        if 'y_cont' not in spec._t.colnames:
            logging.info("I'm adding column 'y_cont'.")
        spec._t['y_cont'] = spec._t['y']
        if 'cont' not in spec._t.colnames:
            logging.info("I'm adding column 'cont'.")
        spec._t['cont'] = at.Column(np.array(None, ndmin=1), dtype=float)
        #spec._t['cont'].unit = spec._t['y'].unit

        sum_sel = len(spec._t)
        for i in range(iter):
            sel = spec._t['y']-spec._t['y_rm']<-kappa*spec._t['dy']
            spec._t['y_em'][sel] = 0
            spec._t['y_abs'][sel] = 1
            x_rm = spec._t['x'][~sel]
            y_rm = running_mean(spec._t['y'][~sel], h=hwindow)

            if i == iter-1 and sum_sel != np.sum(sel):
                logging.warning("Clipping not converged after %i iterations. "
                                "Try iterating more." % (i+1))
            if sum_sel != np.sum(sel):
                sum_sel = np.sum(sel)
            else:
                logging.info("Clipping converged after %i iterations." % (i+1))
                break

            spec._t['y_rm'] = np.interp(spec._t['x'], x_rm, y_rm)
        x_cont = spec._t['x'][np.where(spec._t['y_em'])]
        y_cont = spec._t['y'][np.where(spec._t['y_em'])]
        spec._t['y_cont'] = np.interp(spec._t['x'], x_rm, y_rm)*spec._t['y'].unit
        if lines is not None:
            lines._t['y_cont'] = np.interp(lines._t['x'], x_rm, y_rm)*lines._t['y'].unit

        spec._gauss_convolve(std=std, input_col='y_cont', output_col='cont')

        return 0



    def lines_find(self, std_start=100.0, std_end=0.0, col='y', kind='min',
                   kappa_peaks=5.0, resol=resol_def, append=True):
        """ @brief Find lines
        @details Create a line list by convolving a spectrum with different
        gaussian profiles and finding the peaks in the convolved spectrum
        @param std_start Start standard deviation of the gaussian (km/s)
        @param std_end End standard deviation of the gaussian (km/s)
        @param col Column to convolve
        @param kind Kind of extrema ('min' or 'max')
        @param kappa_peaks Number of standard deviations
        @param resol Resolution
        @param append Append lines to existing line list
        @return 0
        """

        try:
            std_start = float(std_start)
            std_end = float(std_end)
            kappa_peaks = float(kappa_peaks)
            resol = None if resol in [None, 'None'] else float(resol)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if col not in self.sess.spec.t.colnames:
            logging.error(msg_col_miss(col))
            return 0

        #check, resol = resol_check(self.sess.spec, resol)
        if resol is not None:
            logging.info("I'm adding column 'resol'.")
            self.sess.spec._t['resol'] = resol

        #for i, std in enumerate(log2_range(std_start, std_end, -1)):
        self._peaks_found = False
        for i, std in enumerate(np.arange(std_start, std_end, -5)):
            col_conv = col+'_conv'
            self.gauss_convolve(std=std, input_col=col, output_col=col_conv)
            self.peaks_find(col=col_conv, kind='min', kappa=kappa_peaks,
                            append=append or (i>0 and self._peaks_found))

        return 0


    def lines_update(self):
        """ @brief Update line list
        @details Update line list after the systems have been fitted, copying
        fitting parameters and computing the line FWHM. The recipe only works if
        the systems were extracted from the line list that is to be updated.
        @return 0
        """

        systs = self.sess.systs
        lines = self.sess.lines
        spec = self.sess.spec

        systs_x = np.array([])

        # Create new columns in line list (if needed)
        cols = np.append(['syst_id'],
                         np.append(systs._t.colnames[0:2],
                                   systs._t.colnames[3:9]))
        for c in cols:
            if c not in lines._t.colnames:
                logging.info("I'm adding column '%s'." % c)
                if c in ['syst_id']:
                    lines._t[c] = at.Column(np.array(np.nan, ndmin=1), dtype=int)
                elif c in ['func']:
                    lines._t[c] = at.Column(np.array('None', ndmin=1), dtype='S5')
                elif c in ['series']:
                    lines._t[c] = at.Column(np.array('None', ndmin=1), dtype='S100')
                else:
                    lines._t[c] = at.Column(np.array(np.nan, ndmin=1), dtype=float)
        lines._t['fwhm'] = at.Column(np.array(np.nan, ndmin=1), dtype=float)

        count = 0
        for s in systs._t:
            for trans in trans_parse(s['series']):
                #systs_x = np.append(systs_x,
                #                    to_x(s['z0'], trans).to(xunit_def).value)
                systs_x = to_x(s['z0'], trans).to(xunit_def).value
                diff_x = np.abs(lines.x.to(xunit_def).value - systs_x)
                if diff_x.min() == 0:
                    count += 1
                    for c in cols:
                        if c == 'syst_id':
                            lines._t[c][diff_x.argmin()] = s['id']
                        elif c == 'series':
                            lines._t[c][diff_x.argmin()] = trans
                        else:
                            lines._t[c][diff_x.argmin()] = s[c]

        if count==0:
            logging.warning("I couldn't copy the fitting parameters and "
                            "removed the columns. Please check that the "
                            "systems have been extracted from the line list "
                            "that is to be updated.")
            lines._t.remove_columns(np.append(cols, 'fwhm'))
            return 0

        for l in lines._t:
            if l['series'] != 'None':
                fwhm = self._voigt_fwhm(l)
                xpix = np.median(spec._t['xmax']-spec._t['xmin'])
                l['fwhm'] = fwhm

        return 0

    def nodes_cont(self, delta_x=500, kappa_nodes=5.0,
                   smooth=0):
        """ @brief Continuum from nodes
        @details Estimate a continuum by extracting, cleaning, and interpolating
        nodes from regions not affected by lines
        @param delta_x Size of slices (km/s)
        @param kappa_nodes Number of standard deviation away from the window average
        @param smooth Smoothing of the spline
        @return 0
        """

        try:
            delta_x = float(delta_x)
            kappa_nodes = float(kappa_nodes)
            smooth = float(smooth)
        except:
            logging.error(msg_param_fail)
            return 0

        self.nodes_extract(delta_x, au.km/au.s)
        self.nodes_clean(kappa_nodes)
        self.nodes_interp(smooth)

        return 0
