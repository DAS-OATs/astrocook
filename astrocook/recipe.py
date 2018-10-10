from . import *
from .utils import *
from astropy import units as u
from collections import OrderedDict as od
from copy import deepcopy as dc
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
        
        if name == 'line_cont':
            self.objs = ['line', 'spec']
            self.procs = ['mask', 'smooth_lowess']
            self.modes = ['pass', 'pass']
            self.defaults = {}
            self.omits = {}
            
        if name == 'line_ew':
            self.objs = ['line']
            self.procs = ['ew']
            self.modes = [None]
            self.defaults = {}
            self.omits = {}

        if name == 'line_find':
            self.objs = ['spec', 'spec', 'spec']
            self.procs = ['convolve', 'find_extrema', 'select_extrema']
            self.modes = ['pass', None, None]
            self.defaults = {}
            self.omits = {'l'}

        if name == 'line_resid':
            self.objs = ['syst', 'spec', 'spec', 'line']
            self.procs = ['extract_resid', 'convolve', 'select_extrema',
                          'exts_merge']
            self.modes = ['pass', 'pass', None, None]
            self.defaults = {}
            self.omits = {}

        """
        if name == 'resid_find':
            self.objs = ['model', 'model', 'model']
            self.procs = ['convolve', 'find_extrema', 'select_extrema']
            self.modes = ['pass', None, None]
            self.defaults = {}
            self.omits = {'l'}
        """    
        if name == 'spec_cont':
            self.objs = ['spec']
            self.procs = ['convolve']
            self.modes = ['pass']
            self.defaults = {'gauss_sigma': 1000}
            self.omits = {}
                           
        if name == 'syst_find':
            self.objs = ['line', 'line', 'line']
            self.procs = ['create_z', 'match_z', 'map_z']
            self.modes = [None, None, None]
            self.defaults = {}
            self.omits = {}

        if name == 'syst_fit':
            self.objs = ['syst']
            self.procs = ['fit']
            self.modes = [None]
            self.defaults = {}
            self.omits = {'s'}

            
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
        out = self.execute(**kwargs)
        spec = out.spec
        self.line = Line(acs=out, x=spec._exts_sel['X'],
                         y=spec._exts_sel['Y'], xmin=spec._exts_sel['XMIN'],
                         xmax=spec._exts_sel['XMAX'], dy=spec._exts_sel['DY'])
        self.__dict__.update(self.line.__dict__)        
        return out

    def line_resid(self, **kwargs):
        out = self.execute(**kwargs)

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
    
    
