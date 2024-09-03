from .cookbook_edit_old import CookbookEditOld
from .functions import *
from .message import *
from .vars import *

from astropy import units as au
import logging

class CookbookEdit(CookbookEditOld):
    """ Cookbook of utilities for editing data
    """

    def __init__(self):
        super(CookbookEdit, self).__init__()


    def modify_columns(self):
        """@brief Modify columns 🚧
        @details 🚧
        @url other_cb.html#modify-columns
        """

        return 0


    def import_systs(self):
        """@brief Import system list 🚧
        @details 🚧
        @url other_cb.html#import-system-list
        """

        return 0


    def extract(self):
        """@brief Extract 🚧
        @details 🚧
        @url other_cb.html#extract
        """

        return 0


    def mask(self, shift=0, tell=True, sky=True, cond=''):
        """@brief Mask
        @details Mask telluric lines, sky lines, or spectral regions defined by
        a condition
        @url other_cb.html#mask
        @param shift Shift to the barycentric frame (km/s)
        @param tell Mask telluric lines
        @param sky Mask sky lines
        @param cond Condition
        @return 0
        """

        try:
            shift = float(shift)
            tell = str(tell) == 'True'
            sky = str(sky) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if tell: self.telluric_mask(shift=shift)
        if sky: self.sky_mask(shift=shift)
        if cond!='':
            self.mask_cond(cond=cond, new_sess=False)
            self.sess.spec._t['mask'] = np.logical_not(self.sess.spec._t['mask'])
            mask = np.where(np.array(self.sess.spec._t['mask']==0, dtype=bool))
            self.sess.spec._t['y'][mask] = np.nan

        return 0
