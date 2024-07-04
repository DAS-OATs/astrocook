from .cookbook_absorbers_old import CookbookAbsorbersOld
from .vars import resol_def
from .line_list import LineList

import logging
import numpy as np
from scipy.signal import find_peaks

class CookbookAbsorbers(CookbookAbsorbersOld):

    def __init__(self):
        super(CookbookAbsorbers, self).__init__()


    def find_lines(self, kind='abs', prominence=None, append=True):
        """ @brief Find lines
        @details Find absorption or emission lines, based on their prominence.
        @url absorbers_cb.html#find-lines
        @param kind Kind
        @param prominence Prominence
        @param append Append to existing line list
        @return 0
        """

        try:
            kind = str(kind)
            prominence = None if prominence in [None, 'None'] \
                else float(prominence)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if kind not in ['abs', 'em']:
            logging.error("`kind` should be `abs` or `em`. Aborting.")
            return 0

        spec = self.sess.spec
        fact = -1 if kind=='abs' else 1

        ynorm = fact*(spec._t['y'])
        if prominence is None: prominence = 5*(spec._t['dy'])

        peaks, properties = find_peaks(ynorm, prominence=prominence)
        lines = LineList(row=spec._t[peaks], source='y', kind=kind,
                         xunit=spec._xunit, yunit=spec._yunit, meta=spec._meta)
        lines.append_replace(append, self.sess)

        return 0


    def model_lya(self):
        """@brief Model Ly-a forest 🚧
        @details 🚧
        @url absorbers_cb.html#model-ly-a-forest
        """

        return 0


    def model_metals(self):
        """@brief Model metals 🚧
        @details 🚧
        @url absorbers_cb.html#model-metals
        """

        return 0


    def identify_unknown(self):
        """@brief Identify_unknown_lines 🚧
        @details 🚧
        @url absorbers_cb.html#identify-unknown-lines
        """

        return 0


    def check_systs(self):
        """@brief Check system list 🚧
        @details 🚧
        @url absorbers_cb.html#check-system-list
        """

        return 0
