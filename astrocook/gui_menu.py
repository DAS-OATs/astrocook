from . import *
from .vars import *
from .gui_dialog import *
#from .session import Session
#from .model import Model
from .model_list import ModelList
from astropy.io import ascii
import datetime
import wx

prefix = "GUI:"

class GUIMenu(object):

    def __init__(self,
                 gui):
        self._gui = gui

    def bar(self):
        bar = wx.MenuBar()
        self._file = GUIMenuFile(self._gui)
        self._edit = GUIMenuEdit(self._gui)
        self._view = GUIMenuView(self._gui)
        self._snacks = GUIMenuSnacks(self._gui)
        self._meals = GUIMenuMeals(self._gui)
        self._cook = GUIMenuCook(self._gui)
        bar.Append(self._file._menu, "File")
        bar.Append(self._edit._menu, "Edit")
        bar.Append(self._view._menu, "View")
        #bar.Append(snacks._menu, "Snacks")
        bar.Append(self._snacks._menu, "Ingredients")
        #bar.Append(meals._menu, "Meals")
        bar.Append(self._meals._menu, "Recipes")
        #bar.Append(cook._menu, "Cook")
        return bar

    def _item(self, menu, id, title, event):
        item = wx.MenuItem(menu, id, title)
        self._gui._panel_sess.Bind(wx.EVT_MENU, event, item)
        menu.Append(item)

    def _item_graph(self, menu, id, title, key):
        item = wx.MenuItem(menu, id, title)
        self._gui._panel_sess.Bind(
            wx.EVT_MENU, lambda e: self._on_graph(e, key, item), item)
        menu.Append(item)

    def _item_method(self, menu, id, cond, title, targ, obj=None):
        item = wx.MenuItem(menu, id, title+'...')
        item.Enable(cond)
        self._gui._panel_sess.Bind(
            wx.EVT_MENU,
            lambda e: self._on_dialog(e, title, targ, obj), item)
        menu.Append(item)

    def _on_dialog(self, event, title, attr, obj=None):
        if isinstance(attr, list):
            dlg = GUIDialogMethods(self._gui, title, attr, obj)
        else:
            dlg = GUIDialogMethod(self._gui, title, attr, obj)

    def _on_graph(self, event, key, item):
        sel = self._gui._graph_main._sel
        if key in sel:
            sel.remove(key)
            #item.Check(False)
        else:
            sel.append(key)
            #item.Check(True)
        self._gui._graph_main._refresh(self._gui._sess_items)

class GUIMenuCook(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=6000,
                 **kwargs):
        super(GUIMenuCook, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Cook menu here
        self._item(self._menu, start_id+1, "Fit CIV forest in all QSOs...",
                   self._on_civ_full)
        self._item(self._menu, start_id+2, "Test...", self._on_test)


    def _on_civ_full(self, event):
        targ_list = ascii.read('/data/cupani/CIV/targets_3.csv')

        for l in targ_list:
            time_start = datetime.datetime.now()
            t = l['name']
            zem = l['z']
            xmin = l['lambdamin']
            xmax = l['lambdamax']
            self._gui._panel_sess._on_open('/data/cupani/CIV/reduced/'+t\
                                           +'.fits')
            sess_start = self._gui._sess_sel
            if sess_start.spec.meta['object'] == 'J2123-0050':
                sess = sess_start.extract_region(xmin=xmin, xmax=xmax)
            else:
                sess = sess_start

            sess.convolve_gauss(std=10)
            sess.find_peaks(kappa=3.0)
            sess.lines._t.remove_rows(sess.lines.y == 0)
            if np.mean(sess.spec._t['y'])<1 and np.std(sess.spec._t['y'])<1:
                sess.spec._t['cont'] = [1] * len(sess.spec._t)*sess.spec.y.unit
            if 'cont' not in sess.spec._t.colnames:
                sess.extract_nodes(delta_x=1000)
                sess.interp_nodes()

            sess_reg = sess.extract_region(xmin=xmin, xmax=xmax)
            self._gui._panel_sess._on_add(sess_reg, open=False)
            sess_center = dc(sess_reg)
            sess_center.add_syst_from_lines(z_end=20, maxfev=10)#series='unknown')

            sess_reg.lines.t['x'] = (1+sess_center.systs.t['z'])\
                                    *xem_d['Ly_a'].to(sess_reg.spec.x.unit)
            sess_reg.lines.t['logN'] = sess_center.systs.t['logN']

            #sess_reg.add_syst_from_lines(series='SiII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='SiIV', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            sess_reg.add_syst_from_lines(series='CIV', logN=None, b=20.0,
                                         dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='FeII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='MgII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)

            sess_reg.add_syst_from_resids(chi2r_thres=2.0, logN=13.0, b=10.0,
                                          maxfev=10)

            sess_reg.compl_syst(n=10)#, z_start=2.128, z_end=2.1372)
            sess_reg.add_syst_slide(col='deabs')#, z_start=1.6, z_end=1.61)
            sess_reg.merge_syst()
            self._gui._graph_main._refresh(self._gui._sess_items)
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_'
                          +datetime.date.today().isoformat()+'.xxx')
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_latest.xxx')
            time_end = datetime.datetime.now()
            print("%s; computation time: %s" \
                  % (datetime.datetime.now(), time_end-time_start))

    def _on_test(self, event):
        targ_list = ascii.read('/data/cupani/CIV/targets_prova.csv')

        for l in targ_list:
            time_start = datetime.datetime.now()
            t = l['name']
            zem = l['z']
            xmin = l['lambdamin']
            xmax = l['lambdamax']
            #xmin = 381
            #xmax = 382.5
            #xmin = 592
            #xmax = 598
            #xmin = 625
            #xmax = 630
            #xmin = 480
            #xmax = 490
            self._gui._panel_sess._on_open('/data/cupani/CIV/reduced/'+t\
                                           +'.fits')
            sess_start = self._gui._sess_sel
            if sess_start.spec.meta['object'] == 'J2123-0050':
                sess = sess_start.extract_region(xmin=xmin, xmax=xmax)
            else:
                sess = sess_start

            sess.convolve_gauss(std=10)
            sess.find_peaks(kappa=3.0)
            sess.lines._t.remove_rows(sess.lines.y == 0)
            if np.mean(sess.spec._t['y'])<1 and np.std(sess.spec._t['y'])<1:
                sess.spec._t['cont'] = [1] * len(sess.spec._t)*sess.spec.y.unit
            if 'cont' not in sess.spec._t.colnames:
                sess.extract_nodes(delta_x=1000)
                sess.interp_nodes()

            sess_reg = sess.extract_region(xmin=xmin, xmax=xmax)
            self._gui._panel_sess._on_add(sess_reg, open=False)
            #"""
            sess_center = dc(sess_reg)
            sess_center.add_syst_from_lines(z_end=20, maxfev=10)#series='unknown')
            sess_reg.lines.t['x'] = (1+sess_center.systs.t['z'])\
                                    *xem_d['Ly_a'].to(sess_reg.spec.x.unit)
            sess_reg.lines.t['logN'] = sess_center.systs.t['logN']
            #"""
            #sess_reg.add_syst_from_lines(series='SiII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='SiIV', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            sess_reg.add_syst_from_lines(series='CIV', logN=None, b=20.0,
                                         dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='FeII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='MgII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)

            sess_reg.add_syst_from_resids(chi2r_thres=1.0, logN=13.0, b=2.0,
                                          maxfev=10)
            """
            sess_reg.compl_syst(n=10)#, z_start=2.128, z_end=2.1372)
            sess_reg.add_syst_slide(col='deabs')#, z_start=1.6, z_end=1.61)
            sess_reg.merge_syst()
            """
            self._gui._graph_main._refresh(self._gui._sess_items)
            """
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_'
                          +datetime.date.today().isoformat()+'.xxx')
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_latest.xxx')
            """
            time_end = datetime.datetime.now()
            print("%s; computation time: %s" \
                  % (datetime.datetime.now(), time_end-time_start))


class GUIMenuEdit(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=2000,
                 **kwargs):
        super(GUIMenuEdit, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Edit menu here
        self._item_method(self._menu, start_id+301, True, "Extract region",
                          'extract_region')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+311, True, "Convert x axis",
                          'convert_x')
        self._item_method(self._menu, start_id+312, True, "Convert y axis",
                          'convert_y')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+321, True, "Shift to rest frame",
                          'shift_to_rf')
        self._item_method(self._menu, start_id+322, True, "Shift from rest frame",
                          'shift_from_rf')


class GUIMenuFile(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=1000,
                 **kwargs):
        super(GUIMenuFile, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()
        self._start_id = start_id

        # Add items to File menu here
        self._item(self._menu, start_id+1, "Open...\tCtrl+O",
                   lambda e: self._on_open(e, **kwargs))
        self._menu.AppendSeparator()
        self._item_method(
            self._menu, start_id+101,
            len(self._gui._sess_item_sel)>1,
            "Combine sessions...", 'combine', obj=self._gui._panel_sess)
        self._item(self._menu, start_id+102, "Save...\tCtrl+S",
                   lambda e: self._on_save(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+401, "Quit\tCtrl+Q",
                   self._gui._panel_sess._on_close)

    def _on_combine(self, event):
        self._gui._panel_sess._combine()

    def _on_open(self, event, path='.'):
        """ Behaviour for Session > Open """

        wildcard = "Astrocook sessions (*.acs)|*.acs|" \
                   "FITS files (*.fits)|*.fits"
        with wx.FileDialog(self._gui._panel_sess, "Open file", path,
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            path = fileDialog.GetPath()
            name = path.split('/')[-1].split('.')[0]
            print(prefix, "I'm loading session %s..." % path)
            sess = Session(path=path, name=name)
            self._gui._panel_sess._on_add(sess, open=True)

    def _on_save(self, event, path='.'):
        """ Behaviour for Session > Save """

        name = self._gui._sess_sel.name
        with wx.FileDialog(self._gui._panel_sess, "Save session", path, name,
                           wildcard="Astrocook session (*.acs)|*.acs",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            path = fileDialog.GetPath()
            dir = fileDialog.GetDirectory()
            print(prefix, "I'm saving session %s..." % path)
            self._gui._sess_sel.save(path)
            """
            try:
                acs = self
                self.IO.acs_write(self.acs, name, dir)

            except IOError:
                wx.LogError("Cannot save session '%s'." % name)
            """


class GUIMenuMeals(GUIMenu):
    def __init__(self,
                 gui,
                 start_id=5000,
                 **kwargs):
        super(GUIMenuMeals, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Meals menu here
        self._item_method(self._menu, start_id, True, "Find lines",
                          ['convolve_gauss', 'find_peaks'])
        self._item_method(self._menu, start_id+1, True, "Guess continuum",
                          ['convolve_gauss', 'find_peaks', 'extract_nodes',
                           'interp_nodes'])
        self._item_method(self._menu, start_id+2, True, "Fit systems",
                          ['add_syst_from_lines', 'add_syst_from_resids',
                           'add_syst_slide', 'compl_syst'])

class GUIMenuSnacks(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=4000,
                 **kwargs):
        super(GUIMenuSnacks, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Snacks menu here
        self._item_method(self._menu, start_id+101, True, "Convolve with gaussian",
                          'convolve_gauss')
        self._item_method(self._menu, start_id+102, True, "Find peaks", 'find_peaks')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+201, True, "Rebin spectrum",
                          'rebin')
        self._item_method(self._menu, start_id+202, True, "Extract nodes",
                          'extract_nodes')
        self._item_method(self._menu, start_id+203, True, "Interpolate nodes",
                          'interp_nodes')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+301, True, "Add and fit a system",
                          'add_syst')
        self._item_method(self._menu, start_id+302, True, "Add and fit systems from "
                          "line list", 'add_syst_from_lines')
        self._item_method(self._menu, start_id+303, True, "Add and fit systems from "
                          "residuals", 'add_syst_from_resids')
        self._item_method(self._menu, start_id+304, True, "Test and fit systems "
                          "by sliding along spectrum", 'add_syst_slide')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+401, True, "Simulate a system",
                          'simul_syst')
        self._item_method(self._menu, start_id+402, True, "Estimate completeness "
                          "with simulated systems", 'compl_syst')


class GUIMenuView(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=3000,
                 **kwargs):
        super(GUIMenuView, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to View menu here
        self._item(self._menu, start_id+1, "Spectrum table",
                   lambda e: self._on_tab(e, 'spec'))
        self._item(self._menu, start_id+2, "Line table",
                   lambda e: self._on_tab(e, 'lines'))
        self._item(self._menu, start_id+3, "System table",
                   lambda e: self._on_tab(e, 'systs'))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+101, "System detection correctness",
                   lambda e: self._on_ima(e, 'corr'))
        self._item(self._menu, start_id+101, "System detection completeness",
                   lambda e: self._on_ima(e, 'compl'))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+201, "Toggle log x axis", self._on_logx)
        self._item(self._menu, start_id+202, "Toggle log y axis", self._on_logy)
        self._menu.AppendSeparator()
        self._submenu = wx.Menu()
        self._item_graph(self._submenu, start_id+301, "Spectrum",
                         'spec_x_y')
        self._item_graph(self._submenu, start_id+302, "Spectrum error",
                         'spec_x_dy')
        self._item_graph(self._submenu, start_id+303, "Convolved spectrum",
                         'spec_x_conv')
        self._item_graph(self._submenu, start_id+304, "Line list",
                         'lines_x_y')
        self._item_graph(self._submenu, start_id+305, "Masked spectrum",
                         'spec_x_ymask')
        self._item_graph(self._submenu, start_id+306, "Nodes",
                         'spec_nodes_x_y')
        self._item_graph(self._submenu, start_id+307, "Continuum",
                         'spec_x_cont')
        self._item_graph(self._submenu, start_id+308, "Model",
                         'spec_x_model')
        self._item_graph(self._submenu, start_id+309, "De-absorbed",
                         'spec_x_deabs')
        self._item_graph(self._submenu, start_id+310, "Spectral format",
                         'spec_form_x')
        self._item_graph(self._submenu, start_id+311, "System list",
                         'systs_z_series')
        self._menu.AppendSubMenu(self._submenu,  "Toggle graph elements")
        self._item(self._menu, start_id+401, "Toggle normalization",
                   self._on_norm)

    def _on_ima(self, event, obj):
        method = '_ima_'+obj
        getattr(self._gui, method)._on_view(event)

    def _on_logx(self, event):
        self._gui._graph_main._logx = ~self._gui._graph_main._logx
        self._gui._graph_main._refresh(self._gui._sess_items)

    def _on_logy(self, event):
        self._gui._graph_main._logy = ~self._gui._graph_main._logy
        self._gui._graph_main._refresh(self._gui._sess_items)

    def _on_norm(self, event):
        self._gui._graph_main._norm = ~self._gui._graph_main._norm
        self._gui._graph_main._refresh(self._gui._sess_items)

    def _on_tab(self, event, obj):
        method = '_tab_'+obj
        getattr(self._gui, method)._on_view(event)
