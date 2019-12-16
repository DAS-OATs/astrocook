from .functions import *
from .message import *
from .vars import *
from astropy import units as au
from scipy.interpolate import UnivariateSpline as uspline

class CookbookContinuum(object):
    """ Cookbook of utilities for continuum fitting
    """

    def __init__(self):
        pass

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

        spec = self.sess.spec
        if col not in spec.t.colnames:
            logging.error("The spectrum has not a column named '%s'. Please "\
                          "pick another one." % col)
            return 0

        peaks = spec._peaks_find(col, kind, kappa)
        source = [col]*len(peaks.t)
        from .line_list import LineList
        lines = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.y, peaks.dy,
                         source, spec._xunit, spec._yunit, meta=spec._meta)

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

        check, resol = resol_check(self.sess.spec, resol)
        if check:
            logging.info("I'm adding column 'resol'.")
            self.sess.spec._t['resol'] = resol

        #for i, std in enumerate(log2_range(std_start, std_end, -1)):
        for i, std in enumerate(np.arange(std_start, std_end, -5)):
            col_conv = col+'_conv'
            self.gauss_convolve(std=std, input_col=col, output_col=col+'_conv')
            self.peaks_find(col=col_conv, kind='min', kappa=kappa_peaks,
                            append=append or i>0)

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
