from .graph import Graph
from .gui_dialog import GUIDialogMini
from .syst_list import SystList
from .vars import graph_sel
from matplotlib import pyplot as plt
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
from matplotlib.figure import Figure
import numpy as np
import wx

class GUIGraphMain(wx.Frame):
    """ Class for the GUI spectrum graph frame """

    def __init__(self,
                 gui,
                 #sess,
                 title="Spectrum",
                 size_x=wx.DisplaySize()[0]*0.96,
                 size_y=wx.DisplaySize()[1]*0.5,
                 main=True,
                 **kwargs):
        """ Constructor """

        self._gui = gui
        self._title = title
        self._size_x = size_x
        self._size_y = size_y
        self._sel = graph_sel

        self._logx = False
        self._logy = False
        self._norm = False
        self._closed = False
        if main:
            self._gui._graph_main = self
        self._init(**kwargs)
        self.SetPosition((wx.DisplaySize()[0]*0.02, wx.DisplaySize()[1]*0.40))

    def _init(self, **kwargs):
        super(GUIGraphMain, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))


        self._panel = wx.Panel(self)
        self._graph = Graph(self._panel, self._gui, self._sel, **kwargs)
        self._textbar = wx.StaticText(self._panel, 1, style=wx.ALIGN_RIGHT)
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        box_toolbar.Add(self._graph._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._textbar, 1, wx.RIGHT)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._graph._canvas, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP)
        self._panel.SetSizer(self._box)
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)

        #self._gui._statusbar = self.CreateStatusBar()
        move_id = self._graph._canvas.mpl_connect('motion_notify_event',
                                                  self._graph._on_move)
        click_id = self._graph._canvas.mpl_connect('button_press_event',
                                                   self._graph._on_click)


    def _refresh(self, sess, **kwargs):
        if self._closed:
            self._init()
        self._graph._refresh(sess, self._logx, self._logy, self._norm, **kwargs)
        self.Show()

    #def _on_line_new(self, event):
    #    print(self._click_xy)

    def _on_syst_new(self, event):
        sess = self._gui._sess_sel
        for s in sess._series_sel.split(','):
            sess.cb.syst_new(series=s, z=self._graph._cursor._z)
        self._gui._refresh(init_cursor=True)

    def _on_close(self, event):
        self._closed = True
        self.Destroy()

class GUIGraphDetail(GUIGraphMain):

    def __init__(self,
                 gui,
                 #sess,
                 title="Spectrum detail",
                 size_x=wx.DisplaySize()[0]*0.4,
                 size_y=wx.DisplaySize()[1]*0.4,
                 **kwargs):
        super(GUIGraphDetail, self).__init__(gui, title, size_x, size_y,
                                             main=False, **kwargs)
        self._norm = True
        self._gui._graph_det = self
        self._graph._legend = False
        self.SetPosition((wx.DisplaySize()[0]*0.58, wx.DisplaySize()[0]*0.02))

    def _define_lim(self, x, t=None, xspan=30, ymargin=0.1):#, norm=False):
        if t == None:
            t = self._gui._sess_sel.spec.t
            t = t[np.logical_and(~np.isnan(t['x']), ~np.isnan(t['y']))]
        w = np.argmin(np.abs(t['x']-x))
        xmin = t['x'][max(w-xspan, 0)]
        xmax = t['x'][min(w+xspan, len(t)-1)]
        ysel = t['y'][np.where(np.logical_and(t['x']>xmin, t['x']<xmax))]
        yspan = ymargin*(np.max(ysel)-np.min(ysel))
        xlim = (xmin, xmax)
        #if norm:
        ylim = (-0.2, 1.2)
        #else:
        #    ylim = (np.min(ysel)-yspan, np.max(ysel)+yspan)

        return xlim, ylim


class GUIGraphHistogram(GUIGraphMain):

    def __init__(self,
                 gui,
                 title="Column histogram",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.4,
                 **kwargs):
        #self._col_values = col_values
        super(GUIGraphHistogram, self).__init__(gui, title, size_x, size_y,
                                                main=False, **kwargs)
        self._gui._graph_hist = self
        self.SetPosition((wx.DisplaySize()[0]*0.48, wx.DisplaySize()[0]*0.02))

    def _init(self, **kwargs):
        super(GUIGraphMain, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))


        self._panel = wx.Panel(self)
        #self._graph = Graph(self._panel, self._gui, self._sel, **kwargs)
        self._fig = Figure()
        self._canvas = FigureCanvasWxAgg(self._panel, -1, self._fig)
        self._toolbar = NavigationToolbar2Wx(self._canvas)
        self._toolbar.Realize()
        self._textbar = wx.StaticText(self._panel, 1, style=wx.ALIGN_RIGHT)
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        #box_toolbar.Add(self._graph._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._textbar, 1, wx.RIGHT)
        self._box = wx.BoxSizer(wx.VERTICAL)
        #self._box.Add(self._graph._canvas, 1, wx.EXPAND)
        self._box.Add(self._canvas, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP)
        self._panel.SetSizer(self._box)
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _refresh(self, sess, **kwargs):
        if self._closed:
            self._init()
        """
        self._graph._ax.clear()
        self._graph._init_ax(111)
        self._graph._ax.hist(col_values)
        self._graph._canvas.draw()
        """
        if hasattr(self, '_ax'):
            self._ax.clear()
        self._ax = self._fig.add_subplot(111)
        self._ax.tick_params(top=True, right=True, direction='in')#self._init_ax(111)
        label = self._gui._col_tab.GetColLabelValue(self._gui._col_sel)
        self._ax.set_xlabel(label.replace('\n',' '))
        self._ax.set_ylabel('Frequency')
        self._ax.hist(self._gui._col_values, bins=20, align='left')
        self._canvas.draw()
        self.Show()
