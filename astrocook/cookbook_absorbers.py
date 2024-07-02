from .vars import resol_def
from .line_list import LineList

import numpy as np
from scipy.signal import find_peaks

class CookbookAbsorbers(object):

    def __init__(self):
        super(CookbookAbsorbers, self).__init__()


    def lines(self, col='y', append=True):
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
        @param col Column to convolve
        @param append Append lines to existing line list
        @return 0
        """

        try:
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if col not in self.sess.spec.t.colnames:
            logging.error(msg_col_miss(col))
            return 0

        #check, resol = resol_check(self.sess.spec, resol)
        """
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
        """
        spec = self.sess.spec
        col = 'y'
        fact = -1
        ynorm = fact*(spec._t['y']/spec._t['cont'])
        prominence = 5*(np.nanmedian(spec._t['dy']/spec._t['cont']))
        print(prominence)
        peaks, properties = find_peaks(ynorm, prominence=prominence)
        lines = LineList(spec._t['x'][peaks],
                         spec._t['xmin'][peaks],
                         spec._t['xmax'][peaks],
                         spec._t['y'][peaks],
                         spec._t['dy'][peaks],
                         [col]*len(peaks),
                         spec._xunit, spec._yunit,
                         meta=spec._meta)
        if append and self.sess.lines is not None \
            and len(self.sess.lines.t) > 0:
            self.sess.lines._append(lines)
            self.sess.lines._clean()
        else:
            self.sess.lines = lines

        return 0
