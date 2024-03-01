from .functions import *
from .functions import _voigt_par_convert, _fadd
from .message import *
from .vars import *
from astropy import table as at
from astropy import units as au
#from scipy.interpolate import UnivariateSpline as uspline
from copy import deepcopy as dc
import cProfile
from scipy.optimize import root_scalar
from scipy.signal import argrelmin
#from scipy.signal import savgol_filter as sgf
from matplotlib import pyplot as plt
import pstats

#from mpl_toolkits import mplot3d


class CookbookContinuum(object):
    """ Cookbook of utilities for continuum fitting
    @details This cookbook contains utilities to detect absorption features in
    the spectrum and remove them, to determine the continuum level before
    absorption.
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

    def lya_corr(self, zem, input_col='y', mode='basic', logN_thres=100,
                 percentile=100):
        """ @brief Correct spectrum for Lyman-alpha opacity
        @details Correct the spectrum flux for the effective Lyman-alpha opacity.

        Depending on `mode`, the correction factor is computed either from a
        simple integration of the H <span style="font-variant:small-caps;">i</span>
        growth function (`basic`), or using the prescriptions by Inoue et al.
        (2014) (`inoue`).

        The `basic` mode require an upper integration limit that can be either
        provided through `logN_thres` or extracted from the distribution of
        H <span style="font-variant:small-caps;">i</span> column densities,
        from the table of fitted absorption systems (if present). In the latter
        case, the column density corresponding to a given `percentile` of the
        distribution is used as the upper integration limit.

        The correction factor is multiplied to column `input_col` of the
        spectrum and the line list (if present), and the result is saved in
        column `[input_col]_taucorr`.
        @param zem Emisson redshift
        @param input_col Column to correct
        @param mode Correction mode ('basic' or 'inoue')
        @param logN_col Threshold for logarithmic column density
        @param percentile Percentile to compute the threshold from the column
        density distribution (only if logN_col is None)
        @return 0
        """

        try:
            zem = float(zem)
            logN_thres = None if logN_thres in [None, 'None'] else float(logN_thres)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        systs = self.sess.systs

        if logN_thres is None:
            try:
                #logN_thres = np.median(self.sess.systs._t['logN'])
                logN_thres = np.percentile(self.sess.systs._t['logN'], percentile)
                logging.info("I estimated the threshold column density from "
                             "the system table: logN_thres = %2.2f." \
                             % logN_thres)
            except:
                logN_thres = 100
                logging.warning("No systems: I couldn't estimate the threshold "
                                "column density. I am using logN_thres = "\
                                "%2.2f." % logN_thres)
        if mode == 'inoue':
            inoue_all = getattr(spec, '_lya_corr_inoue')(zem, input_col,
                                                         apply=False)
            basic_all = getattr(spec, '_lya_corr_basic')(zem, 100, input_col,
                                                         apply=False)
        else:
            inoue_all = np.zeros(len(spec.x))
            basic_all = np.zeros(len(spec.x))
        basic = getattr(spec, '_lya_corr_basic')(zem, logN_thres, input_col,
                                                 apply=False)


        self._lya_corr = 1+(inoue_all-1)*(basic-1)/(basic_all-1)
        """
        plt.plot(spec.x, inoue_all)
        plt.plot(spec.x, basic_all)
        plt.plot(spec.x, basic)
        plt.plot(spec.x, self._lya_corr, linewidth=3)
        """
        #plt.show()
        taucorr = dc(spec._t[input_col])
        taucorr *= self._lya_corr
        spec._t[input_col+'_taucorr'] = taucorr

        w = spec.x.value > (1+zem)*xem_d['Ly_lim'].value
        logging.info("Mean correction bluewards from Lyman limit: %3.2f." \
                     % np.mean(self._lya_corr[w]))
        return 0


    def nodes_clean(self, kappa=5.0):
        """ @brief Clean nodes
        @details Remove outliers from a list of nodes.

        The recipe applies the peak finding algorithm of
        [`peaks_find`](#find-peaks) to the list of nodes, to detect outliers and
        remove them iterately. Outliers are defined as the peaks that are more
        prominent than `kappa` times the value of `dy` with respect to the
        neighbouring nodes.
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

    def nodes_extract(self, delta_x=500, xunit=au.km/au.s, mode='std'):
        """ @brief Extract nodes
        @details Extract nodes from a spectrum, optionally masking lines.

        Depending on `mode`, the recipe can be used either to create a new set
        of nodes to estimate the continuum (`std`) or to place nodes across an
        existing continuum profile (`cont`). In the latter case, the spectrum
        must have a `cont` column beforehand.

        In `std` mode, the recipe splits the spectrum in a number of slices and
        extracts a node for each slice by averaging the values of `x`, `y`, and
        `dy` within it. The node list is formatted like a normal spectrum.

        Slices are defined in either wavelength or velocity space, as
        specified by the chosen `xunit`. They are equally spaced, with their
        size determined by `delta_x`.

        If the spectrum has a `lines_mask` column, only the regions where this
        mask is `0` are considered when etracting the position of the nodes.
        This means that a slice may not be assigned a node if it is completely
        masked.

        If the spectrum has a `deabs` column, containing the de-absorbed flux
        computed by fitting and removing the absorption features, this column
        is used instead of `y` to determine the `y` position of the nodes.

        In `cont` mode, the `cont` column is used instead of `y` to determine
        the `y` position of the nodes, and the `lines_mask` is ignored.
        @param delta_x Size of slices
        @param xunit Unit of wavelength or velocity
        @param mode Mode ('std' for creating nodes from spectrum, 'cont' for
        placing nodes across existing continuum)
        @return 0
        """
        try:
            xunit = au.Unit(xunit)
            delta_x = float(delta_x)*xunit
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.nodes = self.sess.spec._nodes_extract(delta_x, xunit, mode)
        #print(len(self.sess.nodes.t))

        return 0


    def nodes_interp(self, smooth=0):
        """ @brief Interpolate nodes
        @details Interpolate nodes to estimate the emission level.

        The nodes in the current nodes list are interpolated with
        [`scipy.interpolate.UnivariateSpline`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html).
        The interpolating function is sampled at all values of `x` of the
        spectrum and the line list (if present), and is saved as the current
        continuum level in column `cont`.
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
        @details Find the peaks in a spectrum column.

        Peaks are defined as the extrema (minima or maxima) of `y` column that
        are more prominent than `kappa` times the value of `dy` with respect
        to neighbouring bins.

        The recipe finds all the peaks in the spectrum and saves their `x`,
        `xmin`, `xmax`, `y`, and `dy` values in the current line list. If the
        peak is a minimum, `xmin` and `xmax` are determined as the position of
        the closest maxima, and vice-versa if it is a maximum. If `append` is
        `True`, the new peaks are appended to the existing line list (if
        present), otherwise they replace it.

        A column `line_mask` is created in the spectrum (if not present) and
        set to `1` within the `xmin`-`xmax` range of each detected peak.
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
        #print(len(peaks.t))
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


    def peaks_find_new(self, col='y', window=101, kappa=2):
        """ @brief Find peaks
        @details Find the peaks (i.e the regions with prominent absorption
        features) in a spectral column. Peaks are identified via iterative
        kappa-sigma clipping.
        @param col Column where to look for peaks
        @param window Rolling window for clipping
        @param kappa Number of standard deviations to clip peaks
        @return Selection mask for peaks
        """

        try:
            window = int(window)
            kappa = float(kappa)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        if col not in spec.t.colnames:
            logging.error("The spectrum has not a column named '%s'. Please "\
                          "pick another one." % col)
            return 0

        x = spec.t['x']
        y = spec.t['y']
        dy = spec.t['dy']
        peaks = np.full(y.shape, False)
        speaks = -1
        while np.sum(peaks)!=speaks:
            #if i > 0:
            xpeaks = np.array(x[~peaks])
            ypeaks = np.array(y[~peaks])
            speaks = np.sum(peaks)
            #else:
            #    xpeaks = np.array(x)
            #    ypeaks = np.array(y)
            #    speaks = 0
            csy = np.cumsum(
                      np.insert(np.insert(ypeaks, 0, np.full((window-1)//2+1, ypeaks[0])),
                                -1, np.full((window-1)//2, ypeaks[-1])))
            fy = (csy[window:]-csy[:-window]) / float(window)
            iy = np.interp(x, xpeaks, fy)

            peaks = y<iy-kappa*dy
            #print(np.sum(peaks), speaks)
            #if np.sum(peaks)==speaks: break

            #y = y*1e17
            gy = np.gradient(y)

            #ax = plt.axes(projection ="3d")
            #ax.scatter3D(y[~peaks]-fy[~peaks], np.log(dy[~peaks]), gy[~peaks], s=1)
            #ax.scatter3D(y[peaks]-fy[peaks], np.log(dy[peaks]), gy[peaks], s=2)

            plt.plot(spec.x, iy, color='C1')

        #plt.show()

        return peaks


    def peaks_merge(self, peaks, col='y', append=True):
        """ @brief Merge peaks
        @details Merge the peaks into a list of lines.
        @param peaks Peaks
        @param col Column from which peaks were extracted
        @param append Append peaks to existing line list
        @return 0
        """

        try:
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        #if type(peaks) != type(spec):
        #    logging.error("Peaks must be of Spectrum type.")
        #    return 0
        if col not in spec.t.colnames:
            logging.error("The spectrum has not a column named '%s'. Please "\
                          "pick another one." % col)
            return 0

        dist = 1

        merge = dc(peaks)
        #print(merge)
        for i in range(dist, len(peaks)-dist):
            merge[i] = merge[i]+np.logical_and(np.any(peaks[i-dist:i]), np.any(peaks[i+1:i+dist+1]))
        #print(merge)

        id = merge*(range(len(merge))-np.cumsum(merge))
        #print(id)

        sel = np.array([], dtype=int)
        for i in np.unique(id)[1:]:
            add = argrelmin(spec.t['y'][id==i])[0] if len(spec.t['y'][id==i])>1 else spec.t['y'][id==i]
            sel = np.append(sel, np.array(add, dtype=int)+np.array(range(len(merge)))[id==i][0])
        #print(sel)


        peaks_spec = spec._copy(sel)
        if len(peaks_spec.t) > 0:
            source = [col]*len(peaks_spec.t)
            from .line_list import LineList
            lines = LineList(peaks_spec.x, peaks_spec.xmin, peaks_spec.xmax,
                             peaks_spec.y, peaks_spec.dy,
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


    def flux_clip(self, xmin, xmax, hwindow=100, kappa_hi=6, kappa_lo=3,
                  iter=100, fudge=1, std=500, delta_x=5000):
        """ @brief Clip flux in spectrum
        @details Discriminate absorbed spectrum bins by applying a kappa-sigma
        clipping within a running window.

        The recipes computes the mean of `y` in a runniung window across the
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


    def lines_find(self, std_start=100.0, std_end=0.0, col='y', kind='min',
                   kappa_peaks=5.0, resol=resol_def, append=True):
        """ @brief Find lines
        @details Create a line list by convolving a spectrum with different
        gaussian profiles and finding the peaks in the convolution.

        The recipe is an extension of [peaks_find](#find-peaks), using
        convolution to detect peaks of different width. The gaussian profiles
        are applied in velocity space and used to filter out all fluctuations
        below a given scale.

        The recipe proves most effective when a wide range of line widths is
        spanned, from the widest to the narrowest (typically, with gaussian
        standard deviation ranging from 100 to 0 km/s). A growing list of
        detected peaks is created and updated while spanning the range of line
        widths, eliminating duplicated findings.

        The final list of `x`, `xmin`, `xmax`, `y`, and `dy` values for the
        peaks are saved in the current line list. If `append` is `True`, the new
        peaks are appended to the existing line list (if present), otherwise
        they replace it.

        N.B. If the same line is found more than once, using different line
        widths for the gaussian filtering, it is saved just once in the line
        lists, and the value of `xmin` and `xmax` are the positions of the
        closest extrema (the closes maxima if the peak is a minimum, and
        vice-versa) when the line was first detected. This is the reason for
        spanning line widths in decreasing order: if a wide line is detected at
        first with a small value of gaussian standard deviation, it may be
        assigned a `xmin`-`xmax` range much smaller than the region it actually
        covers (because the filter is too fine to smooth out the noise around
        the peak), and the range is preserved even if the line is detected again
        with wider gaussians.
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


    def lines_find_new(self, col='y', window=301, kappa=3.0, resol=resol_def,
                       append=True):
        """ @brief Find lines
        @details Detect absorption lines by finding peaks and extracting their
        minima.
        @param col Column to use
        @param window Rolling window for clipping
        @param kappa Number of standard deviations to clip peaks
        @param resol Resolution
        @param append Append lines to existing line list
        @return 0
        """

        try:
            window = int(window)
            kappa = float(kappa)
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
        peaks = self.peaks_find_new(col, window, kappa)
        self.peaks_merge(peaks, col, append)

        return 0


    def lines_update(self):
        """ @brief Update line list
        @details Update the line list with the parameters obtained by fitting
        the absorption systems.

        This recipe is meant to be run after the absorption systems have been
        fitted with Voigt profiles. It does not work if the list of systems
        comes from a different line lists.

        For each systems, the best-fitting Voigt parameters are propagated to
        the line list and the FWHM of the lines is computed.
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

    def nodes_cont(self, delta_x=500, kappa_nodes=5.0, smooth=0):
        """ @brief Continuum from nodes
        @details Estimate a continuum by extracting, cleaning, and interpolating
        nodes from regions not affected by lines.

        This recipe combines [nodes_extract](#extract-nodes),
        [nodes_clean](#clean-nodes), and [nodes_interp](#interp-nodes) to
        estimate the continuum level once the lines have been detected.
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

    def abs_cont(self, zem, std=1000.0, resol=resol_def, mode='basic',
                 reest_n=4, _refit_n=0, _percentile=100, _print_stats=True):
        """ @brief Continuum from absorbers
        @details Estimate a continuum by iteratively fitting and removing
        absorbers
        @param zem Emisson redshift
        @param std Standard deviation of the gaussian (km/s)
        @param resol Resolution
        @param mode Correction mode ('basic' or 'inoue')
        @param reest_n Number of re-estimation cycles
        @return 0
        """
        profile = cProfile.Profile()
        profile.enable()
        try:
            zem = float(zem)
            std = float(std)
            resol = None if resol in [None, 'None'] else float(resol)
            reest_n = int(reest_n)
            refit_n = int(_refit_n)
            percentile = float(_percentile)
            print_stats = str(_print_stats) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        lines = self.sess.lines
        systs = self.sess.systs

        self.lya_corr(zem, input_col='y', mode=mode, logN_thres=100)
        self.gauss_convolve(std, input_col='y_taucorr', output_col='cont')
        self.lines_find(resol=resol)
        #print('lines', len(self.sess.lines._t))
        self.systs_new_from_lines(refit_n=refit_n, resol=resol)
        #print('systs', len(self.sess.systs._t))
        for i in range(reest_n):
            spec._t['decorr'] = spec._t['deabs']/self._lya_corr
            spec._t['cont%i' % i] = spec._t['cont']
            spec._t['deabs%i' % i] = spec._t['deabs']
            self.lya_corr(zem, input_col='decorr', mode=mode, logN_thres=None,
                          percentile=percentile)
            self.gauss_convolve(std, input_col='decorr_taucorr',
                                output_col='cont')
            self.lines_find(resol=resol, col='deabs')
            #print('lines', len(self.sess.lines._t))
            if i == reest_n-1 or True:
                self.systs_new_from_lines(refit_n=refit_n, resol=resol, append=False)
            else:
                self.systs_new_from_lines(refit_n=0, resol=resol, append=False)
            #print('systs', len(self.sess.systs._t))
            #self.systs_fit(refit_n=1)
        self.nodes_extract(delta_x=1000.0, mode='cont')
        profile.disable()
        ps = pstats.Stats(profile)
        if print_stats:
            ps.sort_stats('cumtime').print_stats()
        return 0
