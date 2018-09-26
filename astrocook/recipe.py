from . import *
from .utils import *
from astropy import units as u
from collections import OrderedDict as od
from copy import deepcopy as dc
import inspect

class Recipe():
    def __init__(self, obj, name):

        self.obj = obj
        #self.spec = obj.spec
        #self.line = obj.line
        #self.syst = obj.syst
        #self.cont = obj.cont
        self.name = name
        self.descr = rec_descr[name]        
        self.params = None
        
        if name == 'line_find':
            self.procs = ['convolve', 'find_extrema', 'select_extrema']
            self.modes = ['pass', None, None]
            self.defaults = {}

        if name == 'line_cont':
            self.procs = ['mask', 'smooth_lowess']
            self.modes = ['pass', 'pass']
            self.defaults = {}

        if name == 'spec_cont':
            self.procs = ['convolve']
            self.modes = ['pass']
            self.defaults = {'gauss_sigma': 1000}
                           
        if name == 'syst_find':
            self.procs = ['create_z', 'match_z']#, 'create_t']
            self.modes = [None, 'pass']#, None]
            self.defaults = {}

        """
        if name == 'syst_def':
            self.params = {'forest': 'Ly',
                           'zem': 0.0,
                           'prox_vel': 0.0,
                           'xmin': 0.0,
                           'xmax': 0.0}
            self.dialog = od([('Ion:', 'forest'),
                              ('Emission redshift:', 'zem'),
                              ('Prox. velocity:', 'prox_vel'),
                              ('Min. wavelength:', 'xmin'),
                              ('Max. wavelength:', 'xmax')])
            self.procs = ['syst_def']
            
        if name == 'syst_find':
            self.params = {'series': 'CIV',
                           'ztol': 3e-4}
            self.dialog = od([('Series:', 'series'),
                              ('Redshift tolerance:', 'ztol')])
            self.procs = ['find']
        """
        
    def execute(self, **kwargs):

        obj = dc(self.obj)
        for p, m in zip(self.procs, self.modes):
            method = getattr(obj, p)
            try:
                param = {k: kwargs[k] for k in kwargs \
                         if k in inspect.getargspec(method)[0][1:]}
                out = method(**param)
            except:
                out = method()
            if m != None:
                obj = out
        return obj

    def line_cont(self, **kwargs):

        out = self.execute(**kwargs)
        t = out.t
        self.cont = Cont(out, x=t['X'], y=t['Y'], dy=t['DY'])
        self.__dict__.update(self.cont.__dict__)        
        return out

    def line_find(self, **kwargs):

        out = self.execute(**kwargs)
        sel = out._exts_sel
        self.line = Line(self.obj, x=sel['X'], y=sel['Y'], xmin=sel['XMIN'],
                    xmax=sel['XMAX'], dy=sel['DY'])
        self.__dict__.update(self.line.__dict__)        
        return out

    def spec_cont(self, **kwargs):

        out = self.execute(**kwargs)
        t = self.obj.t
        self.cont = Cont(self.obj, x=t['X'], y=t['Y'], dy=t['DY'])
        self.__dict__.update(self.cont.__dict__)        

    def syst_find(self, **kwargs):

        out = self.execute(**kwargs)
        print out
        self.syst = System(line=self.obj, z=out)
        return out

    
