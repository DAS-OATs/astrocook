from . import version
from .session import Session
from collections import OrderedDict
from copy import deepcopy as dc
import inspect
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
import numpy as np
from sphinx.util import docstrings as ds
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix

prefix = "GUI:"

class GUI(object):
    """ Class for the GUI. """

    def __init__(self):
        """ Constructor """

        print("┌───────────────────┐")
        print("│ ASTROCOOK 🔭  v%3s │" % version)
        print("└───────────────────┘")
        print("AC: Welcome! Try Session > Open…")
        self._sess_list = []#Session()
        self._sess_sel = None
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
                 gui,
                 title,
                 targ,
                 attr):

        self._gui = gui
        self._gui._dlg_method = self
        self._targ = targ
        self._attr = attr
        self._obj = getattr(self._gui._sess_sel, self._targ)
        self._method = getattr(self._obj, self._attr)
        super(GUIDialogMethod, self).__init__(parent=None, title=title)

        self._get_params()
        self._get_doc()
        panel = wx.Panel(self)
        box = wx.BoxSizer(wx.VERTICAL)
        core = wx.BoxSizer(wx.VERTICAL)

        # Description
        hdr = wx.BoxSizer(wx.HORIZONTAL)
        st = wx.StaticText(panel, 1, label=self._brief)
        hdr.Add(st, 1, 0, border=15)
        core.Add(hdr, flag=wx.BOTTOM, border=15)

        # Parameters
        static = wx.StaticBox(panel, label="Parameters")
        sizer = wx.StaticBoxSizer(static, wx.VERTICAL)
        fgs = wx.FlexGridSizer(len(self._params), 2, 5, 5)
        fgs_add = []
        self._ctrl = []
        for p, d in zip(self._params, self._doc):
            stat = wx.StaticText(panel, -1, label=d+':')
            ctrl = wx.TextCtrl(panel, -1, value=str(self._params[p]))
            fgs_add.append((stat))
            fgs_add.append((ctrl, 1, wx.EXPAND))
            self._ctrl.append(ctrl)
        fgs.AddMany(fgs_add)
        sizer.Add(fgs, proportion=1, flag=wx.ALL|wx.EXPAND, border=5)
        #core.Add(fgs, proportion=1, flag=wx.ALL|wx.EXPAND, border=5)
        core.Add(sizer)
        panel.SetSizer(core)

        # Buttons
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        cancel_button = wx.Button(self, label='Cancel')
        cancel_button.Bind(wx.EVT_BUTTON, self._on_cancel)
        run_button = wx.Button(self, label='Run')
        run_button.Bind(wx.EVT_BUTTON, self._on_run)
        run_button.SetDefault()
        buttons.Add(cancel_button, 0, wx.RIGHT, border=5)
        buttons.Add(run_button, 0)

        box.Add(panel, 0, wx.EXPAND|wx.ALL, border=10)
        box.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        box.SetSizeHints(self)

        self.SetSizer(box)
        self.Centre()
        self.Show()

    def _get_doc(self):
        full = inspect.getdoc(self._method)
        split = full.split('@')
        self._brief = [s[6:-1] for s in split if s[0:5]=='brief'][0]
        self._doc = [s[6:-1].split(' ', 1)[1] for s in split if s[0:5]=='param']

    def _get_params(self):
        keys = inspect.getargspec(self._method)[0][1:]
        defs = inspect.getargspec(self._method)[-1]
        if defs == None:
            defs = []
        values = np.append(['']*(len(keys)-len(defs)), defs)
        self._params = OrderedDict(zip(keys, values))

    def _on_cancel(self, e):
        self.Close()

    def _on_run(self, e):
        for p, c in zip(self._params, self._ctrl):
            pmod = c.GetValue()
            self._params[p] = pmod
        out = self._method(**self._params)
        if out is not None:
            if out is 0:
                self._gui._panel_sess._refresh()
            else:
                new_sess = dc(self._gui._sess_sel)
                setattr(new_sess, self._targ, out)
                self._gui._panel_sess._on_add(e, new_sess, open=False)
            self.Close()

class GUIMenu(object):

    def __init__(self,
                 gui):
        self._gui = gui

    def bar(self):
        bar = wx.MenuBar()
        session = GUIMenuSession(self._gui)
        spectrum = GUIMenuSpectrum(self._gui)
        bar.Append(session._menu, "Session")
        bar.Append(spectrum._menu, "Spectrum")
        return bar

    def _item(self, menu, id, label, event):
        item = wx.MenuItem(menu, id, label)
        self._gui._panel_sess.Bind(wx.EVT_MENU, event, item)
        menu.Append(item)

    def _item_method(self, menu, id, title, targ, attr):
        item = wx.MenuItem(menu, id, title)
        self._gui._panel_sess.Bind(wx.EVT_MENU,
                                 lambda e: self._on_dialog(e, title, targ, attr),
                                 item)
        menu.Append(item)

    def _on_dialog(self, event, title, targ, attr):
        dlg = GUIDialogMethod(self._gui, title, targ, attr)

class GUIMenuSession(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=1000,
                 **kwargs):
        super(GUIMenuSession, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to session menu here
        open = self._item(self._menu, start_id, "Open",
                         lambda e: self._on_open(e, **kwargs))
        quit = self._item(self._menu, start_id+1, "Quit", self._on_quit)

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
            source = fileDialog.GetPath()
            #chosen = fileDialog.GetDirectory()
            print(prefix, "Loading session %s..." % source)
            sess = Session(source)
            sess._source = source
            sess._name = source.split('/')[-1][:-5]
            self._gui._panel_sess._on_add(event, sess, open=True)
        print("AC: Now try Spectrum…")

    def _on_quit(self, event):
        print("AC: Bye!")
        self._gui._panel_sess.Close()
        self._gui._graph_spec.Close()
        self._gui._tab_spec.Close()

class GUIMenuSpectrum(GUIMenu):

    def __init__(self,
                 gui,
                 start_id=2000,
                 **kwargs):
        super(GUIMenuSpectrum, self).__init__(gui)
        self._gui = gui
        self._menu = wx.Menu()

        # Add items to Spectrum menu here
        tab = self._item(self._menu, start_id, "View table", self._on_view)
        self._menu.AppendSeparator()
        self._item_method(self._menu, start_id+100, "Convert wavelengths",
                          'spec', 'convert_x')
        self._item_method(self._menu, start_id+101, "Convolve with gaussian",
                          'spec', 'convolve_gauss')
        self._item_method(self._menu, start_id+102, "Extract region",
                          'spec', 'extract_region')

    def _on_view(self, event):
        self._gui._tab_spec._on_view(event)

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
        self._tab.InsertColumn(4, "# lines", width=100)
        self._tab.InsertColumn(5, "# systems", width=100)
        self._tab.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self._on_veto)
        self._tab.Bind(wx.EVT_LIST_END_LABEL_EDIT, self._on_edit)
        self._tab.Bind(wx.EVT_LIST_ITEM_SELECTED, self._on_select)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        panel.SetSizer(self._box)
        self._menu = GUIMenu(self._gui)
        self.SetMenuBar(self._menu.bar())
        self.Show()

    def _on_add(self, event, sess, open=True):
        self._sel = self._tab.GetItemCount()
        self._tab.insert_string_item(self._sel, "%s (%s)"
                                     % (sess._name, str(self._sel)))
        self._gui._sess_list.append(sess)
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        if open:
            self._gui._sess_sel.open()
        x = sess.spec._safe(sess.spec.x)#.value
        #xunit = sess.spec.x.unit
        self._refresh()
        """
        self._tab.SetItem(self._sel, 1, obj)
        self._tab.SetItem(self._sel, 2, "[%3.2f, %3.2f] %s"
                          % (x[0].value, x[-1].value, x.unit))
        self._tab.SetItem(self._sel, 3, str(len(x)))
        self._gui._graph_spec._refresh(self._gui._sess_sel)
        """
    def _on_edit(self, event):
        self._gui._sess_list[self._sel].spec.meta['object'] = event.GetLabel()

    def _on_select(self, event):
        self._sel = event.GetIndex()
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        name = self._tab.GetItem(self._tab.GetFirstSelected(), 0)
        self._refresh()#self._gui._sess_sel)

    def _on_veto(self, event):
        if event.GetColumn() in [0,2,3,4,5]:
            event.Veto()
        else:
            event.Skip()

    def _refresh(self):
        sess = self._gui._sess_sel
        obj = sess.spec.meta['object']
        x = sess.spec._safe(sess.spec.x)
        self._tab.SetItem(self._sel, 1, obj)
        self._tab.SetItem(self._sel, 2, "[%3.2f, %3.2f] %s"
                          % (x[0].value, x[-1].value, x.unit))
        self._tab.SetItem(self._sel, 3, str(len(x)))
        self._gui._graph_spec._refresh(self._gui._sess_sel)

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
        self._gui = gui
        self._gui._graph_spec = self


        panel = wx.Panel(self)
        self._fig = Figure()
        self._ax = self._fig.add_subplot(111)
        #self._fig.tight_layout(rect=[-0.03, 0.02, 1.03, 1])
        self._plot = FigureCanvasWxAgg(panel, -1, self._fig)
        self._toolbar = NavigationToolbar2WxAgg(self._plot)
        self._toolbar.Realize()
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        box_toolbar.Add(self._toolbar, 1, wx.RIGHT, border=5)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._plot, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP, border=5)
        panel.SetSizer(self._box)
        self.Centre()

    def _refresh(self, sess):
        self._ax.clear()
        self._ax.plot(sess.spec.x, sess.spec.y)
        self._plot.draw()
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

        self._gui = gui
        self._attr = attr
        self._panel = wx.Panel(self)
        self._tab = gridlib.Grid(self._panel)
        self._tab.CreateGrid(0, 0)

    def _on_view(self, event):
        data = getattr(self._gui._sess_sel, self._attr)
        try:
            self._tab.DeleteCols(pos=0, numCols=self._tab.GetNumberCols())
            #self._tab.DeleteRows(pos=0, numRows=self._tab.GetNumberRows())
            print(prefix, "Refreshing table...")
        except:
            print(prefix, "Loading table...")
        coln = len(data.t.colnames)
        rown = len(data.t)
        self._tab.AppendCols(coln)
        self._tab.AppendRows(rown)
        for j, r in enumerate(data.t):
            for i, n in enumerate(data.t.colnames):
                if j == 0:
                    self._tab.SetColSize(i, 150)
                    self._tab.SetColLabelValue(i, "%s\n%s" \
                                              % (n, str(data.t[n].unit)))
                self._tab.SetCellValue(j, i, "%3.5f" % r[n])
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        self._panel.SetSizer(self._box)
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

        self._gui = gui
        self._gui._tab_spec = self
