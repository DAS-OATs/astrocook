from . import *
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
        item = wx.MenuItem(menu, id, title+'…')
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
        self._item(self._menu, start_id+1, "My recipe…", self._on_my_recipe)

    def _on_my_recipe(self, event):
        sess = self._gui._sess_sel
        sess.convolve_gauss()
        sess.find_peaks()
        sess.extract_nodes()
        sess.interp_nodes()
        new_sess = sess.extract_region(xmin=400, xmax=420)
        self._gui._panel_sess._on_add(new_sess, open=False)
        new_sess.add_syst_from_lines(resol=140000)
        new_sess.test_fit_slide(logN=11.5)
        new_sess.add_syst(z=1.6967)
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
        self._item(self._menu, start_id+1, "Open…",
                   lambda e: self._on_open(e, **kwargs))
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+101, "Save…",
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
            name = path.split('/')[-1][:-5]
            print(prefix, "I'm loading session %s..." % path)
            sess = Session(path=path, name=path.split('/')[-1][:-5])
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
                          ['add_fit_from_lines', 'test_fit_slide'])

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
        self._item_method(self._menu, start_id+303, "Test and fit systems "
                          "by sliding along spectrum", 'test_fit_slide')


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
                   lambda e: self._on_view(e, 'spec'))
        self._item(self._menu, start_id+2, "Line table",
                   lambda e: self._on_view(e, 'lines'))
        self._item(self._menu, start_id+3, "System table",
                   lambda e: self._on_view(e, 'systs'))
        self._item(self._menu, start_id+4, "Model table",
                   lambda e: self._on_view(e, 'mods'))
        self._menu.AppendSeparator()
        self._menu.AppendSeparator()
        self._item(self._menu, start_id+101, "Toggle log x axis", self._on_logx)
        self._item(self._menu, start_id+102, "Toggle log y axis", self._on_logy)
        self._menu.AppendSeparator()
        self._submenu = wx.Menu()
        self._item_graph(self._submenu, start_id+201, "Spectrum",
                         'spec_x_y')
        self._item_graph(self._submenu, start_id+202, "Spectrum error",
                         'spec_x_dy')
        self._item_graph(self._submenu, start_id+203, "Convolved spectrum",
                         'spec_x_conv')
        self._item_graph(self._submenu, start_id+204, "Line list",
                         'lines_x_y')
        self._item_graph(self._submenu, start_id+205, "Masked spectrum",
                         'spec_x_ymask')
        self._item_graph(self._submenu, start_id+206, "Nodes",
                         'spec_nodes_x_y')
        self._item_graph(self._submenu, start_id+207, "Continuum",
                         'spec_x_cont')
        self._item_graph(self._submenu, start_id+208, "Model",
                         'spec_x_model')
        self._item_graph(self._submenu, start_id+208, "De-absorbed",
                         'spec_x_deabs')
        self._item_graph(self._submenu, start_id+209, "Spectral format",
                         'spec_form_x')
        self._menu.AppendSubMenu(self._submenu,  "Toggle graph elements")
        self._item(self._menu, start_id+210, "Toggle normalization",
                   self._on_norm)

    def _on_logx(self, event):
        self._gui._graph_spec._logx = ~self._gui._graph_spec._logx
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_logy(self, event):
        self._gui._graph_spec._logy = ~self._gui._graph_spec._logy
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_norm(self, event):
        self._gui._graph_spec._norm = ~self._gui._graph_spec._norm
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_view(self, event, obj):
        method = '_tab_'+obj
        getattr(self._gui, method)._on_view(event)
