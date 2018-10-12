from . import *
from .utils import *
from .procedure import *
from astropy import units as u
#from collections import OrderedDict as od
from copy import deepcopy as dc
#from copy import copy
import inspect

rec_descr = {'line_cont': "Estimate continuum by masking lines",
             'line_find': "Find lines",
             'line_ew': "Estimate all equivalent widths",
             'line_resid': "Add lines from model residuals to the selected "\
                           "system",
             'spec_cont': "Estimate continuum by smoothing spectrum",
             'syst_find': "Find systems",             
             }



class Recipe(Procedure):

    def __init__(self, acs=None):
        """ @brief Constructor for an abstract recipe
        
        @param acs Session
        """
        
        super(Recipe, self).__init__(acs)

    def ex(self, **kwargs):

        acs = self.acs
        objs = self.objs
        procs = self.procs
        forw = self.forw
        copy = self.copy
        paste = self.paste
        for o, p, cp, fw, ps in zip(objs, procs, copy, forw, paste):
            obj = getattr(acs, o)
            proc = getattr(obj, p)
            if cp != None:
                bck = dc(getattr(acs, cp)) 
            try:
                param = {k: kwargs[k] for k in kwargs \
                         if k in inspect.getargspec(proc)[0][1:]}
                out = proc(**param)
            except:
                pass
            if ps != None:
                setattr(acs, ps, bck)
            if fw != None:
                setattr(acs, fw, out)

    def get_params(self):
        for o, p in zip(self.objs, self.procs):
            self.obj = o
            self.proc = p
            super(Recipe, self).get_params()
                
class RecLineCont(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor of the recipe that estimates continuum by masking
        lines
        
        @param acs Session
        """
        
        super(RecLineCont, self).__init__(acs)
        self.title = rec_descr['line_cont']
        self.objs = ['line', 'spec', 'cont']
        self.procs = ['mask', 'smooth_lowess', 'spec_new']
        self.copy = ['spec', None, None]
        self.forw = ['spec', 'spec', 'cont']
        self.paste = [None, None, 'spec']
        self.defaults = {}
        self.omits = {}
        

class RecLineFind(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that finds lines
        
        @param acs Session
        """

        super(RecLineFind, self).__init__(acs)
        self.title = rec_descr['line_find']
        self.objs = ['spec', 'spec', 'line']
        self.procs = ['convolve', 'select_extrema', 'exts_new']
        self.copy = ['spec', None, None]
        self.forw = ['spec', None, 'line']
        self.paste = [None, None, 'spec']
        self.defaults = {}
        self.omits = {}


class RecLineResid(Recipe):
    
    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that adds lines from model 
        residual to the selected system
        
        @param acs Session
        """

        super(RecLineResid, self).__init__(acs)
        self.title = rec_descr['line_resid']
        self.objs = ['syst', 'spec', 'spec', 'line']
        self.procs = ['extract_resid', 'convolve', 'select_extrema',
                      'exts_merge']
        self.copy = ['spec', None, None, None]
        self.forw = ['spec', 'spec', None, None]
        self.paste = [None, None, None, 'spec']
        self.defaults = {}
        self.omits = {}

    
class RecSpecCont(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that estimates continuum by 
        smoothing spectrum
        
        @param acs Session
        """

        super(RecSpecCont, self).__init__(acs)
        self.title = rec_descr['spec_cont']
        self.objs = ['spec', 'cont']
        self.procs = ['convolve', 'spec_new']
        self.copy = ['spec', None]
        self.forw = ['spec', 'cont']
        self.paste = [None, 'spec']
        self.defaults = {'gauss_sigma': 1000}
        self.omits = {}

        
class RecSystFind(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that finds systems
        
        @param acs Session
        """

        super(RecSystFind, self).__init__(acs)
        self.title = rec_descr['syst_find']
        self.objs = ['line', 'line', 'line', 'syst']
        self.procs = ['create_z', 'match_z', 'map_z', 'line_new']
        self.copy = [None, None, None, None]
        self.forw = [None, None, None, 'syst']
        self.paste = [None, None, None, None]
        self.defaults = {}
        self.omits = {}
    
