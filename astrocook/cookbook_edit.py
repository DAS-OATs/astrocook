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


    def mask(self):
        """@brief Mask 🚧
        @details 🚧
        @url other_cb.html#mask
        """

        return 0
