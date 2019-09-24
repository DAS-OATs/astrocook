from . import * #version
from .graph import Graph
from .gui_image import *
from .gui_menu import *
from .gui_table import *
import numpy as np
from sphinx.util import docstrings as ds
import wx
import wx.lib.mixins.listctrl as listmix

prefix = "GUI:"

class GUI(object):
    """ Class for the GUI. """

    def __init__(self, path=None):
        """ Constructor """

        print("┌───────────────────┐")
        print("│ ASTROCOOK 🍪  v%3s │" % version)
        print("└───────────────────┘")
        #print("-----------------")
        #print(" ASTROCOOK  v%3s " % version)
        #print("-----------------")
        print("Cupani et al. 2017-2019 * INAF-OATs")
        self._sess_list = []
        self._sess_sel = None
        self._panel_sess = GUIPanelSession(self)
        GUIGraphSpectrum(self)
        GUITableSpectrum(self)
        GUITableLineList(self)
        GUITableSystList(self)
        GUITableModelList(self)
        GUIImageCompleteness(self)
        GUIImageCorrectness(self)
        if path == None:
            print("AC: Welcome! Try Session > Open...")
        else:
            self._panel_sess._on_open(path)


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
        self._gui = gui
        self._gui._panel_sess = self

        # Create table
        panel = wx.Panel(self)
        self._tab = GUIControlList(panel, 0)
        self._tab.InsertColumn(0, "name", width=200)
        self._tab.InsertColumn(1, "object", width=150)
        self._tab.InsertColumn(2, "active range", width=200)
        self._tab.InsertColumn(3, "# rows", width=100)
        self._tab.InsertColumn(4, "# nodes", width=100)
        self._tab.InsertColumn(5, "# lines", width=100)
        self._tab.InsertColumn(6, "# systems", width=100)
        self._tab.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self._on_veto)
        self._tab.Bind(wx.EVT_LIST_END_LABEL_EDIT, self._on_edit)
        self._tab.Bind(wx.EVT_LIST_ITEM_SELECTED, self._on_select)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        panel.SetSizer(self._box)
        self._menu = GUIMenu(self._gui)
        self.SetMenuBar(self._menu.bar())
        self.Show()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _on_add(self, sess, open=True):
        # _sel is the last selection; _items is the list of all selections.
        self._sel = self._tab.GetItemCount()
        self._items = [self._sel]

        self._tab.insert_string_item(self._sel, "%s (%s)"
                                     % (sess.name, str(self._sel)))
        self._gui._sess_list.append(sess)

        # Similarly, _sess_sel contains the last selected session; _sess_items
        # contains all selected sessions
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_items = [self._gui._sess_sel]
        if open:
            self._gui._sess_sel.open()
        x = sess.spec._safe(sess.spec.x)#.value
        self._refresh()
        self._gui._graph_spec._refresh(self._gui._sess_items)
        #print(self._gui._sess_sel.__dict__)

    def _on_edit(self, event):
        self._gui._sess_list[self._sel].spec.meta['object'] = event.GetLabel()

    def _on_open(self, path):
        """ Behaviour for Session > Open """

        #name = path.split('/')[-1][:-5]
        name = path.split('/')[-1].split('.')[0]
        print(prefix, "I'm loading session %s..." % path)
        sess = Session(path=path, name=name)
        self._gui._panel_sess._on_add(sess, open=True)

    def _on_close(self, event):
        print("AC: Bye!")
        self.Destroy()
        """
        self._gui._panel_sess.Close()
        self._gui._graph_spec.Close()
        self._gui._tab_spec.Close()
        self._gui._tab_lines.Close()
        self._gui._tab_systs.Close()
        self._gui._tab_mods.Close()
        """
        exit()

    def _on_select(self, event):
        self._sel = event.GetIndex()
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        item = self._tab.GetFirstSelected()
        self._items = []
        while item != -1:
            self._items.append(item)
            item = self._tab.GetNextSelected(item)
        self._gui._sess_items = [self._gui._sess_list[i] for i in self._items]
        self._refresh()
        self._gui._graph_spec._refresh(self._gui._sess_items)

    def _on_veto(self, event):
        if event.GetColumn() in [0,2,3,4,5]:
            event.Veto()
        else:
            event.Skip()

    def _refresh(self):
        for i, s in zip(self._items, self._gui._sess_items):
            obj = s.spec.meta['object']
            self._tab.SetItem(i, 1, obj)
            x = s.spec._safe(s.spec.x)
            self._tab.SetItem(i, 2, "[%3.2f, %3.2f] %s"
                              % (x[0].value, x[-1].value, x.unit))
            self._tab.SetItem(i, 3, str(len(x)))
            try:
                x = s.nodes._safe(s.nodes.x)
                self._tab.SetItem(i, 4, str(len(x)))
            except:
                pass
            try:
                x = s.lines._safe(s.lines.x)
                self._tab.SetItem(i, 5, str(len(x)))
            except:
                pass
            try:
                x = s.systs.z
                self._tab.SetItem(i, 6, str(len(x)))
            except:
                pass

class GUIGraphSpectrum(wx.Frame):
    """ Class for the GUI spectrum graph frame """

    def __init__(self,
                 gui,
                 #sess,
                 title="Spectrum",
                 size_x=wx.DisplaySize()[0]*0.9,
                 size_y=wx.DisplaySize()[1]*0.5):
        """ Constructor """

        self._gui = gui
        self._title = title
        self._size_x = size_x
        self._size_y = size_y
        self._gui._graph_spec = self
        self._sel = ['spec_x_y', 'lines_x_y', 'spec_x_cont', 'spec_x_model',
                     'systs_z_series']

        self._logx = False
        self._logy = False
        self._norm = False
        self._closed = False
        self._init()

    def _init(self):
        super(GUIGraphSpectrum, self).__init__(parent=None, title=self._title,
                                               size=(self._size_x, self._size_y))

        panel = wx.Panel(self)
        self._graph = Graph(panel, self._gui, self._sel)
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        box_toolbar.Add(self._graph._toolbar, 1, wx.RIGHT)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._graph._canvas, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP)
        panel.SetSizer(self._box)
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)

        self._gui._statusbar = self.CreateStatusBar()
        move_id = self._graph._canvas.mpl_connect('motion_notify_event',
                                                  self._graph._on_move)


    def _refresh(self, sess):
        if self._closed:
            self._init()
        self._graph._refresh(sess, self._logx, self._logy, self._norm)
        self.Show()

    def _on_close(self, event):
        self._closed = True
        self.Destroy()
