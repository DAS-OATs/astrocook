from .graph import Graph
import wx

class GUIGraphMain(wx.Frame):
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
        self._gui._graph_main = self
        self._sel = ['spec_x_y', 'lines_x_y', 'spec_x_cont', 'spec_x_model',
                     'systs_z_series']

        self._logx = False
        self._logy = False
        self._norm = False
        self._closed = False
        self._init()

    def _init(self):
        super(GUIGraphMain, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))

        panel = wx.Panel(self)
        self._graph = Graph(panel, self._gui, self._sel)
        self._gui._textbar = wx.StaticText(panel, 1, style=wx.ALIGN_RIGHT)
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        box_toolbar.Add(self._graph._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._gui._textbar, 1, wx.RIGHT)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._graph._canvas, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP)
        panel.SetSizer(self._box)
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)

        #self._gui._statusbar = self.CreateStatusBar()
        move_id = self._graph._canvas.mpl_connect('motion_notify_event',
                                                  self._graph._on_move)


    def _refresh(self, sess, **kwargs):
        if self._closed:
            self._init()
        self._graph._refresh(sess, self._logx, self._logy, self._norm, **kwargs)
        self.Show()

    def _on_close(self, event):
        self._closed = True
        self.Destroy()

class GUIGraphDetail(GUIGraphMain):

    def __init__(self,
                 gui,
                 #sess,
                 title="Spectrum detail",
                 size_x=wx.DisplaySize()[0]*0.4,
                 size_y=wx.DisplaySize()[1]*0.6):
        super(GUIGraphDetail, self).__init__(gui, title, size_x, size_y)
        self._gui._graph_det = self

    def _define_lim(self, row, span=10):
        xmin = row['x']+span*(row['xmin']-row['x'])
        xmax = row['x']+span*(row['xmax']-row['x'])
        return (xmin, xmax), (xmin, xmax)
