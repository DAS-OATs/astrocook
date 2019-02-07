from . import Session
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
import wx
import wx.lib.mixins.listctrl as listmix

class GUI(object):
    """ Class for the GUI. """

    def __init__(self):
        """ Constructor """

        #self.menu = GUIMenu(self)
        GUISessionPanel(self)
        GUISpectrumGraph(self)

class GUIControlList(wx.ListCtrl, listmix.TextEditMixin):
    """ Class for editable control lists."""

    def __init__(self,
                 parent,
                 ID=wx.ID_ANY,
                 pos=wx.DefaultPosition,
                 size=wx.DefaultSize,
                 style=wx.LC_REPORT):
        """ Constructor """

        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)

    def insert_string_item(self, *args):
        self.InsertItem(*args)
        listmix.TextEditMixin.__init__(self)

class GUIMenu(object):

    def __init__(self,
                 gui):
        self.gui = gui

    def bar(self):
        bar = wx.MenuBar()
        file = GUIMenuFile(self.gui)
        bar.Append(file.menu, "&File")
        return bar

    def item(self, menu, id, label, event):
        item = wx.MenuItem(menu, id, label)
        self.gui.sess_panel.Bind(wx.EVT_MENU, event, item)
        menu.Append(item)

class GUIMenuFile(GUIMenu):

    def __init__(self,
                 gui,
                 **kwargs):
        super(GUIMenuFile, self).__init__(gui)
        self.gui = gui

        self.menu = wx.Menu()

        # Add items to File menu here
        open = self.item(self.menu, 100, "&Open\tCtrl+O",
                         lambda e: self.on_open(e, **kwargs))
        quit = self.item(self.menu, 200, "&Quit\tCtrl+Q", None)

    def on_open(self, event, path='.'):
        """ Behaviour for File > Open """

        wildcard = "Astrocook sessions (*.acs)|*.acs|" \
                   "FITS files (*.fits)|*.fits"
        with wx.FileDialog(self.gui.sess_panel, "Open file", path,
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            name = fileDialog.GetPath()
            #chosen = fileDialog.GetDirectory()
            self.gui.sess = Session(name)
            self.gui.sess_panel.on_open(event, name)

class GUISessionPanel(wx.Frame):
    """ Class for the GUI session panel """

    def __init__(self,
                 gui,
                 title="Sessions",
                 size_x=wx.DisplaySize()[0]*0.6,
                 size_y=wx.DisplaySize()[1]*0.2):
        """ Constructor """

        super(GUISessionPanel, self).__init__(parent=None, title=title,
                                             size=(size_x, size_y))

        # Import GUI
        self.gui = gui
        self.gui.sess_panel = self

        # Create table
        self.list = []
        panel = wx.Panel(self)
        self.tab = GUIControlList(panel, 0)
        self.tab.InsertColumn(0, "source", width=150)
        self.tab.InsertColumn(1, "object", width=150)
        self.tab.InsertColumn(2, "redshift", width=150)
        self.tab.InsertColumn(3, "active range [nm]", width=150)
        self.tab.InsertColumn(4, "# lines", width=150)
        self.tab.InsertColumn(5, "# systems", width=150)
        self.tab.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_select)
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.box.Add(self.tab, 1, wx.EXPAND)
        panel.SetSizer(self.box)
        self.menu = GUIMenu(self.gui)
        self.SetMenuBar(self.menu.bar())
        self.Show()

    def on_open(self, event, name):
        self.sel = self.tab.GetItemCount()
        self.tab.insert_string_item(self.sel, name)
        self.list.append(self.gui.sess)
        self.list[self.sel].open()
        self.gui.spec_graph.refresh()

    def on_select(self, event):
        self.sel = event.GetIndex()
        name = self.tab.GetItem(self.tab.GetFirstSelected(), 0)
        self.gui.sess = self.list[self.sel]
        self.gui.spec_graph.refresh()


class GUISpectrumGraph(wx.Frame):
    """ Class for the GUI spectrum graph frame """

    def __init__(self,
                 gui,
                 title="Spectrum",
                 size_x=wx.DisplaySize()[0]*0.9,
                 size_y=wx.DisplaySize()[1]*0.5):
        """ Constructor """

        super(GUISpectrumGraph, self).__init__(parent=None, title=title,
                                              size=(size_x, size_y))
        self.gui = gui
        self.gui.spec_graph = self

        panel = wx.Panel(self)
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        #self.fig.tight_layout(rect=[-0.03, 0.02, 1.03, 1])
        self.plot = FigureCanvasWxAgg(panel, -1, self.fig)
        self.toolbar = NavigationToolbar2WxAgg(self.plot)
        self.toolbar.Realize()
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        box_toolbar.Add(self.toolbar, 1, wx.RIGHT, border=5)
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.box.Add(self.plot, 1, wx.EXPAND)
        self.box.Add(box_toolbar, 0, wx.TOP, border=5)
        panel.SetSizer(self.box)
        self.Centre()

    def refresh(self):
        self.ax.clear()
        self.ax.plot(self.gui.sess.spec.x, self.gui.sess.spec.y)
        self.plot.draw()
        self.Show()
