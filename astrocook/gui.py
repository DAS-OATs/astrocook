from . import version
from .session import Session
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
import numpy as np
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix

prefix = "GUI:"

class GUI(object):
    """ Class for the GUI. """

    def __init__(self):
        """ Constructor """

        print("┌────────────────────────────┐")
        print("│ Welcome to ASTROCOOK v%3s! │" % version)
        print("└────────────────────────────┘")
        print(prefix, "Loading...")
        self.sess_list = []#Session()
        self.sess_sel = None
        GUIPanelSession(self)
        GUIGraphSpectrum(self)
        GUITableSpectrum(self)

class GUIControlList(wx.ListCtrl, listmix.TextEditMixin):
    """ Class for editable control lists. """

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

class GUIDialogMethod(wx.Dialog):
    def __init__(self,
                 parent=None,
                 title="Dialog"):

        super(GUIDialogMethod, self).__init__(parent=parent, title=title)

class GUIMenu(object):

    def __init__(self,
                 gui):
        self.gui = gui

    def bar(self):
        bar = wx.MenuBar()
        file = GUIMenuFile(self.gui)
        spectrum = GUIMenuSpectrum(self.gui)
        bar.Append(file.menu, "&File")
        bar.Append(spectrum.menu, "&Spectrum")
        return bar

    def item(self, menu, id, label, event):
        item = wx.MenuItem(menu, id, label)
        self.gui.panel_sess.Bind(wx.EVT_MENU, event, item)
        menu.Append(item)

class GUIMenuFile(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=1000,
                 **kwargs):
        super(GUIMenuFile, self).__init__(gui)
        self.gui = gui
        self.menu = wx.Menu()

        # Add items to File menu here
        open = self.item(self.menu, start_id, "&Open\tCtrl+O",
                         lambda e: self.on_open(e, **kwargs))
        quit = self.item(self.menu, start_id+1, "&Quit\tCtrl+Q", self.on_quit)

    def on_open(self, event, path='.'):
        """ Behaviour for File > Open """

        wildcard = "Astrocook sessions (*.acs)|*.acs|" \
                   "FITS files (*.fits)|*.fits"
        with wx.FileDialog(self.gui.panel_sess, "Open file", path,
                           wildcard=wildcard,
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            name = fileDialog.GetPath()
            #chosen = fileDialog.GetDirectory()
            print(prefix, "Loading session %s..." % name)
            sess = Session(name)
            sess.name = name
            self.gui.panel_sess.on_open(event, sess)

    def on_quit(self, event):
        self.gui.panel_sess.Close()
        self.gui.graph_spec.Close()
        self.gui.tab_spec.Close()

class GUIMenuSpectrum(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=2000,
                 **kwargs):
        super(GUIMenuSpectrum, self).__init__(gui)
        self.gui = gui
        self.menu = wx.Menu()

        # Add items to Spectrum menu here
        tab = self.item(self.menu, start_id, "View table", self.on_view)

    def on_view(self, event):
        self.gui.tab_spec.on_view(event)

class GUIPanelSession(wx.Frame):
    """ Class for the GUI session panel """

    def __init__(self,
                 gui,
                 title="Sessions",
                 size_x=wx.DisplaySize()[0]*0.6,
                 size_y=wx.DisplaySize()[1]*0.2):
        """ Constructor """

        super(GUIPanelSession, self).__init__(parent=None, title=title,
                                             size=(size_x, size_y))

        # Import GUI
        self.gui = gui
        self.gui.panel_sess = self

        # Create table
        panel = wx.Panel(self)
        self.tab = GUIControlList(panel, 0)
        self.tab.InsertColumn(0, "source", width=200)
        self.tab.InsertColumn(1, "object", width=150)
        self.tab.InsertColumn(2, "active range", width=200)
        self.tab.InsertColumn(3, "# rows", width=100)
        self.tab.InsertColumn(4, "# lines", width=100)
        self.tab.InsertColumn(5, "# systems", width=100)
        self.tab.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self.on_veto)
        self.tab.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.on_edit)
        self.tab.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_select)
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.box.Add(self.tab, 1, wx.EXPAND)
        panel.SetSizer(self.box)
        self.menu = GUIMenu(self.gui)
        self.SetMenuBar(self.menu.bar())
        self.Show()

    def on_edit(self, event):
        self.gui.sess_list[self.sel].spec.meta['object'] = event.GetLabel()

    def on_open(self, event, sess):
        sel = self.tab.GetItemCount()
        self.tab.insert_string_item(sel, sess.name)
        self.gui.sess_list.append(sess)
        self.gui.sess_sel = self.gui.sess_list[sel]
        self.gui.sess_sel.open()
        x = sess.spec.x
        obj = sess.spec.meta['object']
        self.tab.SetItem(sel, 1, obj)
        self.tab.SetItem(sel, 2, "[%3.2f, %3.2f] %s"
                         % (x[~np.isnan(x)][0], x[~np.isnan(x)][-1], x.unit))
        self.tab.SetItem(sel, 3, str(len(x)))
        self.gui.graph_spec.refresh(self.gui.sess_sel)

    def on_select(self, event):
        sel = event.GetIndex()
        self.gui.sess_sel = self.gui.sess_list[sel]
        name = self.tab.GetItem(self.tab.GetFirstSelected(), 0)
        self.gui.graph_spec.refresh(self.gui.sess_sel)

    def on_veto(self, event):
        if event.GetColumn() in [0,2,3,4,5]:
            event.Veto()
        else:
            event.Skip()

class GUIGraphSpectrum(wx.Frame):
    """ Class for the GUI spectrum graph frame """

    def __init__(self,
                 gui,
                 title="Spectrum",
                 size_x=wx.DisplaySize()[0]*0.9,
                 size_y=wx.DisplaySize()[1]*0.5):
        """ Constructor """

        super(GUIGraphSpectrum, self).__init__(parent=None, title=title,
                                              size=(size_x, size_y))
        self.gui = gui
        self.gui.graph_spec = self


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

    def refresh(self, sess):
        self.ax.clear()
        self.ax.plot(sess.spec.x, sess.spec.y)
        self.plot.draw()
        self.Show()

class GUITable(wx.Frame):
    """ Class for the GUI table frame """

    def __init__(self,
                 gui,
                 attr,
                 title="Table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.9):

        super(GUITable, self).__init__(parent=None, title=title,
                                       size=(size_x, size_y))

        self.gui = gui
        self.attr = attr
        self.panel = wx.Panel(self)
        self.tab = gridlib.Grid(self.panel)
        self.tab.CreateGrid(0, 0)

    def on_view(self, event):
        data = getattr(self.gui.sess_sel, self.attr)
        try:
            self.tab.DeleteCols(pos=0, numCols=self.tab.GetNumberCols())
            #self.tab.DeleteRows(pos=0, numRows=self.tab.GetNumberRows())
            print(prefix, "Refreshing table...")
        except:
            print(prefix, "Loading table...")
        coln = len(data.t.colnames)
        rown = len(data.t)
        self.tab.AppendCols(coln)
        self.tab.AppendRows(rown)
        for j, r in enumerate(data.t):
            for i, n in enumerate(data.t.colnames):
                if j == 0:
                    self.tab.SetColSize(i, 150)
                    self.tab.SetColLabelValue(i, "%s\n%s" \
                                              % (n, str(data.t[n].unit)))
                self.tab.SetCellValue(j, i, "%3.5f" % r[n])
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.box.Add(self.tab, 1, wx.EXPAND)
        self.panel.SetSizer(self.box)
        self.Centre()
        self.Show()

class GUITableSpectrum(GUITable):
    """ Class for the GUI spectrum table frame """

    def __init__(self,
                 gui,
                 title="Spectrum table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.9):

        super(GUITableSpectrum, self).__init__(gui, 'spec', title, size_x,
                                               size_y)

        self.gui = gui
        self.gui.tab_spec = self
