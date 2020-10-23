from .message import *
import logging
import numpy as np

class CookbookTemplates(object):
    """ Cookbook of spectral templates
    """

    def __init__(self):
        super(CookbookTemplates, self).__init__()

    def bb(self, temp=6000, scale=1.0):
        """ @brief Blackbody template
        @details Create a blackbody template of a given temperature
        @param temp Temperature (K)
        @param scale Scale factor
        @return 0
        """

        try:
            temp = float(temp)
            scale = float(scale)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._template_bb(temp, scale)
        return 0

    def pl(self, ampl=1.0, x_ref=None, index=1.0):
        """ @brief Power-law template
        @details Create a power-law template of a given amplitude, wavelength
        reference and index
        @param ampl Amplitude
        @param x_ref Wavelength reference (nm)
        @param index Index
        @return 0
        """

        try:
            ampl = float(ampl)
            x_ref = np.mean(self.sess.spec.x) if x_ref in [None, 'None'] \
                else float(x_ref)
            index = float(index)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._template_pl(ampl, x_ref, index)
        return 0
