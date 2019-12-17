from .vars import *
from .cookbook_absorbers import CookbookAbsorbers
from .cookbook_continuum import CookbookContinuum
from .cookbook_general import CookbookGeneral
from .cookbook_sandbox import CookbookSandbox
from .format import Format
from .spectrum import Spectrum
from .syst_list import SystList
from .syst_model import SystModel
from astropy import constants as ac
from astropy import units as au
from copy import deepcopy as dc
import datetime
import logging
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

class Cookbook(CookbookGeneral,
               CookbookContinuum,
               CookbookAbsorbers,
               CookbookSandbox):
    """ Main cookbook, combining specific cookbooks.

    Each cookbook should link to methods of classes containing the actual
    algorithms. All method parameters must have default values."""

    def __init__(self,
                 sess=None):
        super(Cookbook, self).__init__()
        self.sess = sess

    def _refresh(self, sess):
        self.sess = sess
