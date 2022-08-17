#import locale
#locale.setlocale(locale.LC_ALL, 'en_US')
#import wx
#wx.Locale(wx.LANGUAGE_DEFAULT)
version = "1.5.1"
current_year = 2022

from .session import Session
from .syst_list import SystList

import logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(module)s: %(message)s [deprecated]")

#from .frame import Frame
#from .spectrum import Spectrum

#from .format import Format
#from .session import Session


"""
from .list import List  # Deprecated
from .list_reader import ListReader  # Deprecated
from .list_syst import ListSyst  # Deprecated
from .spec_1d import Spec1D
from .spec_1d_cont import Spec1DCont
from .spec_1d_reader import Spec1DReader
from .model import Model
from .line import Line
from .cont import Cont
from .syst import Syst
from .system import System
#from .abs import Abs
from .voigt import Voigt
from .io import IO
from .recipe import Recipe
from .workflow import Workflow
from .plot import Plot
from .app import MainFrame
"""
