from .procedure import *
from .recipe import *
import sys

wkf_descr = {'forest_add_all': "Add new lines of all systems to forest",
             'line_resid_all': "Add residuals of all systems to line list",
             'syst_all_fit': "Fit all systems",
             'syst_all_model': "Model all systems",
             'syst_fit_add': "Fit selected forest system adding lines if "\
                             "needed",
             'syst_fit_add_all': "Fit all forest systems adding lines if "\
                                 "needed",
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

            
    def ex_iter_op(self, op, iter_start=0, iter_end=-1, reverse=False,
                   done=False):

        self.obj = op.obj
        self.proc = op.proc
        obj, proc = self.get_obj_proc()

        if iter_end == -1:
            iter_end = len(obj.t)
        ran = range(iter_start, iter_end)
        if reverse == True:
            ran = reversed(ran)

        for i in ran:
            try:
                if obj.t[i]['DONE'] == False:
                    sys.stdout.write("%i %f " % (i+1, obj.t[i]['Z']))
                    proc(s=obj.t[i], done=done)
                    print obj.t[i]['CHI2R']
            except:
                pass

    def ex_iter_rec(self, rec, iter_start=0, iter_end=-1, reverse=False,
                    done=False):

        root_obj = getattr(self.acs, self.obj)
        root_obj_static = dc(root_obj)
        print root_obj_static.t[0:50]
        if iter_end == -1:
            iter_end = len(root_obj.t)
        ran = range(iter_start, iter_end)
        if reverse == True:
            ran = reversed(ran)

        for i in ran:
            try:
                print i+1, root_obj_static.t[i]['Z']
                root_obj.group(s=root_obj_static.t[i])
                root_obj.chunk()
                rec.ex()
            except:
                pass

                    
    def ex_iter_wkf(self, wkf, iter_start=0, iter_end=-1, reverse=False,
                    done=False):
        
        obj = getattr(self.acs, self.obj)
        if iter_end == -1:
            iter_end = len(obj.t)
        ran = range(iter_start, iter_end)
        if reverse == True:
            ran = reversed(ran)
            
        for i in ran:
            od = False
            for io, o in enumerate(wkf.ops):
                if io == len(wkf.ops)-1:
                    od = done
                if obj.t[i]['DONE'] == False:
                    op = o(acs=self.acs)
                    #print i, od, o, obj.t[i]['Z'], obj.t[i]['DONE']
                    if io == 0:
                        print i+1, obj.t[i]['Z']
                    op.ex(s=obj.t[i], done=od)
                    #op.ex(s=t[i])
                
                
    def get_params(self):
        for o, m, s in zip(self.ops, self.modes, self.setups):
            self.defaults = s
            self.omits = []
            if m == None:
                op = o(acs=self.acs)
                op.get_params()
            else:
                super(Workflow, self).get_params(self, m)


class WkfForestAddAll(Workflow):    

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to add residuals to line list 
        for all systems
        """

        super(WkfForestAddAll, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['forest_add_all']
        self.ops = [RecForestAdd]
        self.obj = 'syst'
        self.modes = [self.ex_iter_rec]
        self.setups = [{'done': True}]
        
                
class WkfLineResidAll(Workflow):    

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to add residuals to line list 
        for all systems
        """

        super(WkfLineResidAll, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['line_resid_all']
        self.ops = [RecLineResid]
        self.obj = 'syst'
        self.modes = [self.ex_iter_rec]
        self.setups = [{}]
        
                
            
class WkfSystFitAll(Workflow):    

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to fit all systems
        """

        super(WkfSystFitAll, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_all_fit']
        self.ops = [ProcSystSelFit]
        self.modes = [self.ex_iter_op]
        self.setups = [{'done': True}]
        

class WkfSystModelAll(Workflow):    

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to model all systems
        """

        super(WkfSystModelAll, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_all_model']
        self.ops = [ProcSystSelModel]
        self.modes = [self.ex_iter_op]
        self.setups = [{'done': True}]
        

class WkfSystFitAdd(Workflow):

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to fit all systems, adding 
        components
        """

        super(WkfSystFitAdd, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_fit_add']
        self.ops = [ProcSystSelFit, RecLineResid, RecForestAdd, ProcLineEwAll,
                    ProcSystNAll, ProcSystSelFit]
        self.modes = [None, None, None, None, None, None]
        self.setups = [{}, {}, {}, {}, {}, {}]


class WkfSystFitAddAll(Workflow):

    def __init__(self, acs=None):
        """ @brief Constructor for the workflow to fit all systems, adding 
        components
        """

        super(WkfSystFitAddAll, self).__init__(acs)
        self.acs = acs
        self.title = wkf_descr['syst_fit_add_all']
        self.ops = [WkfSystFitAdd]
        self.obj = 'syst'
        self.modes = [self.ex_iter_wkf]
        self.setups = [{'done': True}]
    
