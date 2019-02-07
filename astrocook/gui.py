from .gui_elem import GUIControlList
import wx

class GUI(object):
    """ Class for the GUI. """

    def __init__(self):
        """ Constructor """

        self.sess = GUISessions(self)
        self.spec_graph = GUISpectrumGraph(self)

class GUISessions(wx.Frame):
    """ Class for the GUI sessions frame """

    def __init__(self,
                 gui,
                 title="Astrocook",
                 size_x=wx.DisplaySize()[0]*0.6,
                 size_y=wx.DisplaySize()[1]*0.2):
        """ Constructor """

        super(GUISessions, self).__init__(parent=None, title=title,
                                          size=(size_x, size_y))
        panel = wx.Panel(self)
        self.tab = GUIControlList(panel, 0)
        self.tab.InsertColumn(0, 'target', width=150)
        self.tab.InsertColumn(1, 'object', width=150)
        self.tab.InsertColumn(2, 'redshift', width=150)
        self.tab.InsertColumn(3, 'active range [nm]', width=150)
        self.tab.InsertColumn(4, '# lines', width=150)
        self.tab.InsertColumn(5, '# systems', width=150)
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.box.Add(wx.StaticText(panel, label="Sessions"))
        self.box.Add(self.tab, 1, wx.EXPAND)
        panel.SetSizer(self.box)
        self.Centre()
        self.Show()


class GUISpectrumGraph(wx.Frame):
    """ Class for the GUI spectrum graph frame """

    def __init__(self,
                 gui,
                 title="Spectrum",
                 size_x=wx.DisplaySize()[0]*0.6,
                 size_y=wx.DisplaySize()[1]*0.4):
        """ Constructor """

        super(GUISpectrumPlot, self).__init__(parent=None, title=title,
                                              size=(size_x, size_y))
        panel = wx.Panel(self)
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.fig.tight_layout(rect=[-0.03, 0.02, 1.03, 1])
        #self.plot = Plot(self.ax)
        self.Centre()
        self.Show()
