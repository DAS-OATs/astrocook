from . import *
from .vars import *
from .gui_dialog import *
#from .session import Session
#from .model import Model
from .model_list import ModelList
import wx

prefix = "GUI:"

class GUIMenu(object):

    def __init__(self,
                 gui):
        self._gui = gui

    def bar(self):
        bar = wx.MenuBar()
        file = GUIMenuFile(self._gui)
        edit = GUIMenuEdit(self._gui)
        view = GUIMenuView(self._gui)
        snacks = GUIMenuSnacks(self._gui)
        meals = GUIMenuMeals(self._gui)
        cook = GUIMenuCook(self._gui)
        bar.Append(file._menu, "File")
        bar.Append(edit._menu, "Edit")
        bar.Append(view._menu, "View")
        bar.Append(snacks._menu, "Snacks")
        bar.Append(meals._menu, "Meals")
        bar.Append(cook._menu, "Cook")
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

    def _item_method(self, menu, id, title, targ):
        item = wx.MenuItem(menu, id, title+'...')
        self._gui._panel_sess.Bind(
            wx.EVT_MENU,
            lambda e: self._on_dialog(e, title, targ), item)
        menu.Append(item)

    def _on_dialog(self, event, title, attr):
        dlg = GUIDialogMethod(self._gui, title, attr)

    def _on_graph(self, event, key, item):
        sel = self._gui._graph_spec._sel
        if key in sel:
            sel.remove(key)
            #item.Check(False)
        else:
            sel.append(key)
            #item.Check(True)
        self._gui._graph_spec._refresh(self._gui._sess_items)

class GUIMenuCook(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=6000,
                 **kwargs):
        super(GUIMenuCook, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Cook menu here
        self._item(self._menu, start_id+1, "Full CIV search...", self._on_full)
        self._item(self._menu, start_id+2, "Test...", self._on_test)

    def _on_full(self, event):
        test_data_zem = {#'J0003-2323': 2.280,
                         'J0100+0211': 1.959,
                         'J0124-3744': 2.190,
                         'J0240-2309': 2.225,
                         'J1124-1705': 2.39723,
                         'J1344-1035': 2.134,
                         'J1451-2329': 2.208,
                         'J1626+6426': 2.320,
                         #'J2123-0050': 2.26902,
                         }

        if self._gui._sess_sel != None:
            targ = [self._gui._sess_sel.spec.meta['object']]
        else:
            targ = [t for t in test_data_zem]
        for t in targ:

            if self._gui._sess_sel == None:
                self._gui._panel_sess._on_open('test_data/'+t+'.fits')
            sess = self._gui._sess_sel

            zem = test_data_zem[sess.spec.meta['object']]
            xmin = xem_d[series_d['Ly'][-1]].value*(1+zem)
            xmax = xem_d[series_d['CIV'][1]].value*(1+zem)

            sess.convolve_gauss()
            sess.find_peaks()

            if 'cont' not in sess.spec._t.colnames:
                sess.extract_nodes(delta_x=800)
                sess.interp_nodes()
            new_sess = sess.extract_region(xmin=xmin, xmax=xmax)
            self._gui._panel_sess._on_add(new_sess, open=False)
            new_sess.add_syst_from_lines(series='CIV')
            new_sess.add_syst_from_resids(chi2r_thres=1.0, maxfev=100)
            new_sess.add_syst_slide(col='deabs')
            new_sess.compl_syst()

            self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_test(self, event):
        xmin = 360
        xmax = 460
        z_start = 1.32
        z_end = 1.959
        logN_start = 12.0
        logN_end = 11.0
        logN_step = -0.5
        b_start = 5
        b_end = 6
        b_step = 2
        test_data_zem = {'J0003-2323': 2.280,
                         'J0100+0211': 1.959,
                         'J0124-3744': 2.190,
                         'J0240-2309': 2.225,
                         'J1124-1705': 2.39723,
                         'J1344-1035': 2.134,
                         'J1451-2329': 2.208,
                         'J1626+6426': 2.320,
                         'J2123-0050': 2.26902,
                         }
        sess = self._gui._sess_sel
        zem = test_data_zem[sess.spec.meta['object']]
        xmin = xem_d[series_d['Ly'][-1]].value*(1+zem)
        xmax = xem_d[series_d['CIV'][1]].value*(1+zem)
        sess.convolve_gauss()
        sess.find_peaks()
        if 'cont' not in sess.spec._t.colnames:
            sess.extract_nodes(delta_x=800)
            sess.interp_nodes()
        new_sess = sess.extract_region(xmin=xmin, xmax=xmax)
        self._gui._panel_sess._on_add(new_sess, open=False)
        """
        new_sess.add_syst_from_lines(series='CIV', z_start=1.71, z_end=1.18,
                                     resol=70000)
        #new_sess.add_syst_from_resids(chi2r_thres=1.0, maxfev=100)
        #new_sess.add_syst_slide(logN_start=12.0, logN_end=11.0, logN_step=-0.1,
        #                        b_start=2, b_end=10, b_step=2,
        #                        maxfev=100)
        self._gui._graph_spec._refresh(self._gui._sess_items)
        """
        """
        new_sess.add_syst_slide(z_start=z_start, z_end=z_end,
                                logN_start=logN_start, logN_end=logN_end,
                                logN_step=logN_step,
                                b_start=b_start, b_end=b_end, b_step=b_step,
                                maxfev=100)
        """
        """
        new_sess.corr_syst(z_start=z_start, z_end=z_end,
                           logN_start=logN_start, logN_end=logN_end,
                           logN_step=logN_step,
                           b_start=b_start, b_end=b_end, b_step=b_step)
        """
        new_sess.compl_syst(n=100, z_start=z_start, z_end=z_end,
                            logN_start=logN_start, logN_end=logN_end,
                            logN_step=logN_step,
                            b_start=b_start, b_end=b_end, b_step=b_step)
        self._gui._graph_spec._refresh(self._gui._sess_items)


class GUIMenuFile(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=1000,
                 **kwargs):
        super(GUIMenuFile, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to File menu here
        self._item(self._menu, start_id+1, "Open...",
                   lambda e: self._on_open(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+101, "Save...",
                   lambda e: self._on_save(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+401, "Quit", self._on_quit)

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
    def _on_quit(self, event):
        print("AC: Bye!")
        self._gui._panel_sess.Close()
        self._gui._graph_spec.Close()
        self._gui._tab_spec.Close()
        self._gui._tab_lines.Close()
        self._gui._tab_systs.Close()
        self._gui._tab_mods.Close()


class GUIMenuEdit(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=2000,
                 **kwargs):
        super(GUIMenuEdit, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Edit menu here
        self._item_method(self._menu, start_id+301, "Extract region",
                          'extract_region')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+311, "Convert x axis",
                          'convert_x')
        self._item_method(self._menu, start_id+312, "Convert y axis",
                          'convert_y')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+321, "Shift to rest frame",
                          'shift_to_rf')
        self._item_method(self._menu, start_id+322, "Shift from rest frame",
                          'shift_from_rf')


class GUIMenuMeals(GUIMenu):
    def __init__(self,
                 gui,
                 start_id=5000,
                 **kwargs):
        super(GUIMenuMeals, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Meals menu here
        self._item_method(self._menu, start_id, "Find lines",
                          ['convolve_gauss', 'find_peaks'])
        self._item_method(self._menu, start_id+1, "Guess continuum",
                          ['convolve_gauss', 'find_peaks', 'extract_nodes',
                           'interp_nodes'])
        self._item_method(self._menu, start_id+2, "Fit systems",
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
        self._item_method(self._menu, start_id+101, "Convolve with gaussian",
                          'convolve_gauss')
        self._item_method(self._menu, start_id+102, "Find peaks", 'find_peaks')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+201, "Extract nodes",
                          'extract_nodes')
        self._item_method(self._menu, start_id+202, "Interpolate nodes",
                          'interp_nodes')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+301, "Add and fit a system",
                          'add_syst')
        self._item_method(self._menu, start_id+302, "Add and fit systems from "
                          "line list", 'add_syst_from_lines')
        self._item_method(self._menu, start_id+303, "Add and fit systems from "
                          "residuals", 'add_syst_from_resids')
        self._item_method(self._menu, start_id+304, "Test and fit systems "
                          "by sliding along spectrum", 'add_syst_slide')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+401, "Simulate a system",
                          'simul_syst')
        self._item_method(self._menu, start_id+402, "Estimate completeness "
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
        self._gui._graph_spec._logx = ~self._gui._graph_spec._logx
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_logy(self, event):
        self._gui._graph_spec._logy = ~self._gui._graph_spec._logy
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_norm(self, event):
        self._gui._graph_spec._norm = ~self._gui._graph_spec._norm
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_tab(self, event, obj):
        method = '_tab_'+obj
        getattr(self._gui, method)._on_view(event)
