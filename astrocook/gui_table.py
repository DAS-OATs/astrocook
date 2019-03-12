from collections import OrderedDict
import pprint
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix

prefix = "GUI:"

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
            print(prefix, "I'm updating table...")
        except:
            print(prefix, "I'm loading table...")
        coln = len(data.t.colnames)
        rown = len(data.t)-self._tab.GetNumberRows()
        self._tab.AppendCols(coln)
        self._tab.AppendRows(rown)
        for j, r in enumerate(data.t):
            for i, n in enumerate(data.t.colnames):
                if j == 0:
                    self._tab.SetColSize(i, 150)
                    self._tab.SetColLabelValue(i, "%s\n%s" \
                                              % (n, str(data.t[n].unit)))
                try:
                    self._tab.SetCellValue(j, i, "%3.5f" % r[n])
                except:
                    if type(r[n]) == str:
                        self._tab.SetCellValue(j, i, r[n])
                    if type(r[n]) == OrderedDict:
                        self._tab.SetCellValue(j, i, pprint.pformat(r[n]))
                    if type(r[n]) == dict:
                        self._tab.SetCellValue(j, i, pprint.pformat(r[n]))
        self._tab.AutoSizeColumns(True)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        self._panel.SetSizer(self._box)
        self.Centre()
        self.Show()

class GUITableLineList(GUITable):
    """ Class for the GUI line list """

    def __init__(self,
                 gui,
                 title="Line table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.9):

        super(GUITableLineList, self).__init__(gui, 'lines', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_lines = self

class GUITableModelList(GUITable):
    """ Class for the GUI model list """

    def __init__(self,
                 gui,
                 title="Model table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.9):

        super(GUITableModelList, self).__init__(gui, 'mods', title,
                                                size_x, size_y)

        self._gui = gui
        self._gui._tab_mods = self

class GUITableSpectrum(GUITable):
    """ Class for the GUI spectrum table """

    def __init__(self,
                 gui,
                 title="Spectrum table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.9):

        super(GUITableSpectrum, self).__init__(gui, 'spec', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_spec = self

class GUITableSystList(GUITable):
    """ Class for the GUI system list """

    def __init__(self,
                 gui,
                 title="System table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.9):

        super(GUITableSystList, self).__init__(gui, 'systs', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_systs = self
