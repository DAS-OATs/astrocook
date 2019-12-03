from .functions import log2_range
from .message import *
from astropy import units as au
from scipy.interpolate import UnivariateSpline as uspline

class CookbookContinuum(object):
    """ Cookbook of utilities for continuum fitting
    """

    def __init__(self):
        pass

    def nodes_clean(self, kappa=5.0):
        """ @brief Clean nodes
        @details Clean the list of nodes from outliers.
        @param window Number of nodes in the window (must be odd)
        @param kappa Number of standard deviation away from the window average
        @return 0
        """
        try:
            kappa = float(kappa)
        except:
            logging.error(msg_param_fail)
            return 0

        """
        if window%2 == 0:
            window += 1
            logging.warning("window was not odd. I incremented it by one.")
        """

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
            append = append or append == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        if col not in spec.t.colnames:
            logging.error("The spectrum has not a column named '%s'. Please "\
                          "pick another one." % col)
            return 0

        peaks = spec._peaks_find(col, kind, kappa)

        from .line_list import LineList
        lines = LineList(peaks.x, peaks.xmin, peaks.xmax, peaks.y, peaks.dy,
                         spec._xunit, spec._yunit, spec._meta)

        if append and self.sess.lines != None:
            self.sess.lines._append(lines)
            self.sess.lines._clean()
        else:
            self.sess.lines = lines

        self.sess.lines_kind = 'peaks'
        spec._lines_mask(self.sess.lines)
        return 0


    def lines_find(self, std_start=100.0, std_end=0.0, kind='min', kappa=5.0):
        """ @brief Find peaks
        @details Find the peaks in a spectrum column. Peaks are the extrema
        (minima or maxima) that are more prominent than a given number of
        standard deviations. They are saved as a list of lines.
        @param std_start Start standard deviation of the gaussian (km/s)
        @param std_end End standard deviation of the gaussian (km/s)
        @param kind Kind of extrema ('min' or 'max')
        @param kappa Number of standard deviations
        @return 0
        """

        try:
            std_start = float(std_start)
            std_end = float(std_end)
            kappa = float(kappa)
        except:
            logging.error(msg_param_fail)
            return 0

        #for i, std in enumerate(log2_range(std_start, std_end, -1)):
        for i, std in enumerate(np.arange(std_start, std_end, -5)):
            self.gauss_convolve(std=std)
            self.peaks_find(kind='min', kappa=kappa, append=i>0)

        return 0
