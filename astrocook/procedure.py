import inspect

proc_descr = {'convolve': "Convolve with a custom profile",
              'extract_forest': "Extract forest",
              'extract_mask': "Extract non masked regions into new spectrum",
              'extract_resid': "Extract residuals into new spectrum",
              'extract_reg': "Extract spectral region",
              'line_ew_all': "Estimate all equivalent widths",
              'syst_sel_fit': "Fit selected system",
              'mask': "Mask lines",
              #'model': "Model selected system",
              'syst_sel_model': "Model selected system",
              'N_all': "Estimate all column densities",
              'select_extrema': "Select the most prominent extrema",
              'smooth_lowess': "Smooth with a LOWESS method",
}

class Procedure(object):

    def __init__(self, acs=None):
        """ @brief Constructor for an abstract procedure

        @ param acs Session
        """

        self.acs = acs
        self.params = {}

    def ex(self, **kwargs):
        obj, proc = self.get_obj_proc()
        try:
            ok
        except:
            if kwargs != {}:
                param = {k: kwargs[k] for k in kwargs \
                         if k in inspect.getargspec(proc)[0][1:]}
                out = proc(**param)
            else:
                out = proc()            
        #except:
        #    raise Exception("Procedure %s failed!" % self.name)

    def get_obj_proc(self):
        obj = getattr(self.acs, self.obj)
        proc = getattr(obj, self.proc)
        return obj, proc
        
    def get_params(self, obj=None, proc=None):
        if obj == None and proc == None:
            obj, proc = self.get_obj_proc()
        try:
            defs = inspect.getargspec(proc)[3]
            if defs != None:
                keys = inspect.getargspec(proc)[0][-len(defs):]
                self.params.update({k: d for (d, k) in zip(defs, keys)})
            for d in self.defaults:
                if d in self.params:
                    self.params[d] = self.defaults[d] 
            for om in self.omits:
                if om in self.params:
                    del self.params[om]
        except:
            raise Exception("Can't get parameters of %s!" % self.name)

class ProcLineEwAll(Procedure):
    def __init__(self, acs=None):
        """ @brief Constructor for the procedure that estimates all equivalent
        widths

        @ param acs Session
        """

        super(ProcLineEwAll, self).__init__(acs)
        self.name = 'line_ew_all'
        self.title = proc_descr[self.name]
        self.obj = 'line'
        self.proc = 'ew_all'
        self.defaults = {}
        self.omits = {}
    
        

class ProcSystSelModel(Procedure):
    
    def __init__(self, acs=None):
        """ @brief Constructor for the procedure that models selected system

        @ param acs Session
        """

        super(ProcSystSelModel, self).__init__(acs)
        self.name = 'syst_sel_model'
        self.title = proc_descr[self.name]
        self.obj = 'syst'
        self.proc = 'model'
        self.defaults = {}
        self.omits = {}


class ProcSystSelFit(Procedure):
    
    def __init__(self, acs=None):
        """ @brief Constructor for the procedure that fits selected system

        @ param acs Session
        """

        super(ProcSystSelFit, self).__init__(acs)
        self.name = 'syst_sel_fit'
        self.title = proc_descr[self.name]
        self.obj = 'syst'
        self.proc = 'fit'
        self.defaults = {}
        self.omits = {}
        
        
