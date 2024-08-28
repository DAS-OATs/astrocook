from .cookbook_absorbers_old import CookbookAbsorbersOld
from .functions import resol_check, to_z, trans_parse
from .line_list import LineList
from .vars import resol_def, xem_d

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


    def model_metals(self, series, zem, no_ly=True, use_lines=False):
        """@brief Model metals
        @details Model metal absorbers, based on transition and emission
        redshift.
        @url absorbers_cb.html#model-metals
        @param series Transitions
        @param zem Emission redshift
        @param no_ly Exclude Lyman forest
        @param use_lines Use line list to define absorbers
        @url absorbers_cb.html#model-metals
        """

        try:
            zem = float(zem)
            no_ly = str(no_ly) == 'True'
            use_lines = str(use_lines) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if no_ly:
            trans = trans_parse(series)
            z_start = to_z(xem_d['Ly_a']*(1+zem), trans[0])
        else:
            z_start = 0
        z_end = zem

        check, resol = resol_check(self.sess.spec)
        if resol is None:
            self.sess.spec._resol_est(3, True)
            resol = self.sess.spec._t['resol'][0]

        col = 'deabs' if 'deabs' in self.sess.spec._t.colnames else 'y'

        if use_lines and self.sess.lines is None:
            logging.warning("I didn't find a line list. Ignoring it.")
            use_lines = False

        if use_lines:
            self.systs_new_from_lines(series=series, z_start=z_start,
                                      z_end=z_end, refit_n=0, append=True)
        else:
            self.systs_new_from_like(series=series, col=col, z_start=z_start,
                                     z_end=z_end, modul=10, sigma=3,
                                     distance=3, resol=resol, append=True)
        self.systs_fit(refit_n=0)
        return 0


    def identify_unknown(self):
        """@brief Identify unknown lines 🚧
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
