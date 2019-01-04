from .procedure import *
from .recipe import *

wkf_descr = {'syst_all_fit': "Fit all systems",
             'syst_all_model': "Model all systems",
             'syst_fit_add': "Fit selected system adding lines if needed",
}

class Workflow(Procedure):

    def __init__(self, acs=None):
        """ @brief Constructor for an abstract workflow
        """

        super(Workflow, self).__init__(acs)
        
    def ex(self, **kwargs):
        for o, m, s in zip(self.ops, self.modes, self.setups):
            op = o(acs=self.acs)
            if m == None:
                op.ex()
            else:
                param = {k: kwargs[k] for k in kwargs \
                         if k in inspect.getargspec(m)[0][1:]}
                m(op, **param)

            
    def ex_iter_table(self, op, iter_start=0, iter_end=-1):

        self.obj = op.obj
        self.proc = op.proc
        obj, proc = self.get_obj_proc()
        for r in obj.t[iter_start:iter_end]:
            proc(s=r)

                
    def get_params(self):
        for m, s in zip(self.modes, self.setups):
            if m != None:
                self.defaults = s
                self.omits = []
                super(Workflow, self).get_params(self, m)
            

            
class WkfSystAllFit(Workflow):    

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to fit all systems
        """

        super(WkfSystAllFit, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_all_fit']
        self.ops = [ProcSystSelFit]
        self.modes = [self.ex_iter_table]
        self.setups = [{}]
        

class WkfSystAllModel(Workflow):    

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to model all systems
        """

        super(WkfSystAllModel, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_all_model']
        self.ops = [ProcSystSelModel]
        self.modes = [self.ex_iter_table]
        self.setups = [{}]
        

class WkfSystFitAdd(Workflow):

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to fit all systems, adding 
        components
        """

        super(WkfSystFitAdd, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_fit_add']
        self.ops = [ProcSystSelFit, RecLineResid, ProcLineEwAll, ProcSystSelFit]
        self.modes = [None, None, None, None]
        self.setups = [{}, {}, {}, {}]
        
    
