from .gui_dialog import *
from .session import Session
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
        self._gui._panel_sess.Bind(wx.EVT_MENU,
                                   lambda e: self._on_graph(e, key, item), item)
        menu.Append(item)

    def _item_method(self, menu, id, title, source, targ, attr):
        item = wx.MenuItem(menu, id, title+'…')
        self._gui._panel_sess.Bind(
            wx.EVT_MENU,
            lambda e: self._on_dialog(e, title, source, targ, attr), item)
        menu.Append(item)

    def _on_dialog(self, event, title, source, targ, attr):
        dlg = GUIDialogMethod(self._gui, title, source, targ, attr)

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
        pass

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
            self._gui._panel_sess._on_add(event, sess, open=True)

    def _on_quit(self, event):
        print("AC: Bye!")
        self._gui._panel_sess.Close()
        self._gui._graph_spec.Close()
        self._gui._tab_spec.Close()
        self._gui._tab_lines.Close()


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
                          None, None, 'extract_region')
        self._item_method(self._menu, start_id+302, "Convert x axis",
                          None, None, 'convert_x')
        self._item_method(self._menu, start_id+303, "Convert y axis",
                          None, None, 'convert_y')


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
                          ['spec', 'spec'], [None, 'lines'],
                          ['convolve_gauss', 'find_peaks'])
        self._item_method(self._menu, start_id+1, "Guess continuum",
                          ['spec', 'spec', 'spec', 'spec'],
                          [None, 'lines', 'nodes', None],
                          ['convolve_gauss', 'find_peaks', 'extract_nodes',
                           'interp_nodes'])

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
                          'spec', None, 'convolve_gauss')
        self._item_method(self._menu, start_id+102, "Find peaks",
                          'spec', 'lines', 'find_peaks')
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+201, "Extract nodes",
                          'spec', 'nodes', 'extract_nodes')
        self._item_method(self._menu, start_id+202, "Interpolate nodes",
                          'spec', None, 'interp_nodes')

    def _on_view(self, event):
        self._gui._tab_spec._on_view(event)

class GUIMenuView(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=3000,
                 **kwargs):
        super(GUIMenuView, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to View menu here
        self._item(self._menu, start_id+1, "View spectrum table",
                   lambda e: self._on_view(e, 'spec'))
        self._item(self._menu, start_id+2, "View line table",
                   lambda e: self._on_view(e, 'lines'))
        self._item(self._menu, start_id+3, "View system table",
                   lambda e: self._on_view(e, 'systems'))
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
        self._item_graph(self._submenu, start_id+208, "Spectral format",
                         'spec_x_model')
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
