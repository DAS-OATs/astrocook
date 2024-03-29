from .vars import *
from .cookbook_absorbers import CookbookAbsorbers
from .cookbook_continuum import CookbookContinuum
from .cookbook_edit import CookbookEdit
from .cookbook_flux import CookbookFlux
from .cookbook_general import CookbookGeneral
from .cookbook_graph import CookbookGraph
from .cookbook_sandbox import CookbookSandbox
from .cookbook_synthetic import CookbookSynthetic
from .cookbook_templates import CookbookTemplates
from .cookbook_view import CookbookView
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

class Cookbook(CookbookEdit,
               CookbookView,
               CookbookGeneral,
               CookbookContinuum,
               CookbookAbsorbers,
               CookbookFlux,
               CookbookGraph,
               CookbookSynthetic,
               CookbookTemplates):
    """ Main cookbook, combining specific cookbooks.

    Each cookbook should link to methods of classes containing the actual
    algorithms. All method parameters must have default values."""

    def __init__(self,
                 sess=None):
        super(Cookbook, self).__init__()
        self.sess = sess
        self._tag = "cb"


    def _refresh(self, sess):
        self.sess = sess
