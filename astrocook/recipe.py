from . import *
from .utils import *
from astropy import units as u
from collections import OrderedDict as od
from copy import deepcopy as dc
from copy import copy
import inspect

class Recipe():
    def __init__(self, acs, name):

        self.acs = acs
        self.spec = acs.spec
        try:
            self.line = acs.line
        except:
            self.line = None
        try:
            self.syst = acs.syst
        except:
            self.line = None
        try:
            self.cont = acs.cont
        except:
            self.line = None

        self.name = name
        self.descr = rec_descr[name]        
        self.params = None
        
        if name == 'line_resid':
            self.objs = ['syst', 'spec', 'spec', 'line']
            self.procs = ['extract_resid', 'convolve', 'select_extrema',
                          'exts_merge']
            self.modes = ['pass', 'pass', None, None]
            self.defaults = {}
            self.omits = {}
                           

    def ex(self, **kwargs):

        acs = self.acs
        objs = self.objs
        nexts = self.objs[:-1]+[self.objs[-1]]
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
                ok
            except:
                param = {k: kwargs[k] for k in kwargs \
                         if k in inspect.getargspec(proc)[0][1:]}
                out = proc(**param)
            #except:
            #    pass
            if ps != None:
                setattr(acs, ps, bck)
            if fw != None:
                setattr(acs, fw, out)

                        
                           
    def execute(self, **kwargs):

        #acs = dc(self.acs)
        acs = self.acs
        pf = False
        for o, p, m in zip(self.objs, self.procs, self.modes):
            if pf:
                setattr(acs, o, out)
            obj = getattr(acs, o)
            method = getattr(obj, p)
            try:
            #    ok
            #except:
                param = {k: kwargs[k] for k in kwargs \
                         if k in inspect.getargspec(method)[0][1:]}
                out = method(**param)
            except:
                out = method()
            if m != None:
                pf = True
                setattr(acs, o, out)
            else:
                pf = False
        return acs

    def line_cont(self, **kwargs):

        out = self.execute(**kwargs)
        spec = out.spec
        self.cont = Cont(spec=spec, x=spec.t['X'], y=spec.t['Y'],
                         dy=spec.t['DY'])
        self.__dict__.update(self.cont.__dict__)        
        return out

    def line_ew(self, **kwargs):

        for l in self.acs.line.t:
            kwargs['l'] = l
            out = self.execute(**kwargs)
        return out
    
    def line_find(self, **kwargs):
        acs = self.execute(**kwargs)
        self.line = acs.line
        """
        spec = out.spec
        self.line = Line(acs=out, x=spec._exts_sel['X'],
                         y=spec._exts_sel['Y'], xmin=spec._exts_sel['XMIN'],
                         xmax=spec._exts_sel['XMAX'], dy=spec._exts_sel['DY'])
        self.__dict__.update(self.line.__dict__)        
        """
        return acs

    def line_resid(self, **kwargs):
        self.execute(**kwargs)

    def spec_cont(self, **kwargs):

        out = self.execute(**kwargs)
        spec = out.spec
        self.cont = Cont(spec=spec, x=spec.t['X'], y=spec.t['Y'],
                         dy=spec.t['DY'])
        self.__dict__.update(self.cont.__dict__)        

    def syst_find(self, **kwargs):
        out = self.execute(**kwargs)
        line = out.line
        self.syst = System(acs=out, series=kwargs['series'], z=line._z_match)
        print self.syst._map
        return out

    def syst_fit(self, **kwargs):
        out = self.execute(**kwargs)
        syst = out.syst
        #self.syst = System(acs=out, series=kwargs['series'], z=line._z_match)
        return out
    
    
class RecipeLineCont(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor of the recipe that estimates continuum by masking
        lines
        
        @param acs Session
        """

        self.acs = acs
        self.title = rec_descr['line_cont']
        self.objs = ['line', 'spec', 'cont']
        self.procs = ['mask', 'smooth_lowess', 'spec_new']
        self.copy = ['spec', None, None]
        self.forw = ['spec', 'spec', 'cont']
        self.paste = [None, None, 'spec']
        self.defaults = {}
        self.omits = {}
        

class RecipeLineFind(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that finds lines
        
        @param acs Session
        """

        self.acs = acs
        self.title = rec_descr['line_find']
        self.objs = ['spec', 'spec', 'line']
        self.procs = ['convolve', 'select_extrema', 'exts_new']
        self.copy = ['spec', None, None]
        self.forw = ['spec', None, 'line']
        self.paste = [None, None, 'spec']
        self.defaults = {}
        self.omits = {'l'}


class RecipeLineResid(Recipe):
    
    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that adds lines from model 
        residual to the selected system
        
        @param acs Session
        """

        self.acs = acs
        self.title = rec_descr['line_resid']
        self.objs = ['syst', 'spec', 'spec']#, 'line']
        self.procs = ['extract_resid', 'convolve', 'select_extrema']#,'exts_merge']
        self.copy = ['spec', None, None]#, None]
        self.forw = ['spec', 'spec', None]#, None]
        self.paste = [None, None, None]#, 'spec']
        self.defaults = {}
        self.omits = {}

    
class RecipeSpecCont(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that estimates continuum by 
        smoothing spectrum
        
        @param acs Session
        """

        self.acs = acs
        self.title = rec_descr['spec_cont']
        self.objs = ['spec', 'cont']
        self.procs = ['convolve', 'spec_new']
        self.copy = ['spec', None]
        self.forw = ['spec', 'cont']
        self.paste = [None, 'spec']
        self.defaults = {'gauss_sigma': 1000}
        self.omits = {}

        
class RecipeSystFind(Recipe):

    def __init__(self, acs=None):
        """ @brief Constructor for the recipe that finds systems
        
        @param acs Session
        """

        self.acs = acs
        self.title = rec_descr['syst_find']
        self.objs = ['line', 'line', 'line', 'syst']
        self.procs = ['create_z', 'match_z', 'map_z', 'line_new']
        self.copy = [None, None, None, None]
        self.forw = [None, None, None, 'syst']
        self.paste = [None, None, None, None]
        self.defaults = {}
        self.omits = {}
    
