from .vars import *
from collections import OrderedDict
import logging
import matplotlib.pyplot as plt
import numpy as np
import pprint
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix

class GUITable(wx.Frame):
    """ Class for the GUI table frame """

    def __init__(self,
                 gui,
                 attr,
                 title="Table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.2):

        self._gui = gui
        self._attr = attr
        self._title = title
        self._size_x = size_x
        self._size_y = size_y
        super(GUITable, self).__init__(parent=None, title=self._title,
                                       size=(self._size_x, self._size_y))

    def _fill(self):
        for j, r in enumerate(self._data.t):
            for i, n in enumerate(self._data.t.colnames):
                if j == 0:
                    self._tab.SetColSize(i, 150)
                    self._tab.SetColLabelValue(i, "%s\n%s" \
                                              % (n, str(self._data.t[n].unit)))
                if type(r[n]) == np.int64:
                    self._tab.SetCellValue(j, i, "%4i" % r[n])
                elif type(r[n]) == str or type(r[n]) == np.str_:
                    self._tab.SetCellValue(j, i, r[n])
                elif type(r[n]) == OrderedDict:
                    self._tab.SetCellValue(j, i, pprint.pformat(r[n]))
                elif type(r[n]) == dict:
                    self._tab.SetCellValue(j, i, pprint.pformat(r[n]))
                else:
                    self._tab.SetCellValue(j, i, "%3.5f" % r[n])
        self._tab.AutoSizeColumns(True)


    def _init(self, from_scratch=True):
        if from_scratch:
            super(GUITable, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))

            self._panel = wx.Panel(self)
            self._tab = gridlib.Grid(self._panel)
            self._tab.CreateGrid(0, 0)
            self.SetPosition((0, wx.DisplaySize()[1]*0.5))
        else:
            if self._tab.GetNumberRows() != 0:
                self._tab.DeleteRows(0, self._tab.GetNumberRows())

        coln = len(self._data.t.colnames)
        rown = len(self._data.t)-self._tab.GetNumberRows()
        self._tab.AppendCols(coln)
        self._tab.AppendRows(rown)


    def _on_close(self, event):
        self.Destroy()


    def _on_detail(self, event):
        if event.GetRow() == -1: return
        if not hasattr(self._gui, '_graph_det'):
            from .gui_graph import GUIGraphDetail
            GUIGraphDetail(self._gui, init_ax=False)
        elif len(self._gui._graph_det._graph._fig.axes) > 1:
            self._gui._graph_det._graph._fig.clear()
        size_x = wx.DisplaySize()[0]*0.4
        size_y = wx.DisplaySize()[1]*0.4
        self._gui._graph_det.SetSize(wx.Size(size_x, size_y))
        self._gui._graph_det._graph._init_ax(111)
        #row = self._data.t[self._gui._tab_popup._event.GetRow()]
        row = self._data.t[event.GetRow()]
        self._gui._sess_sel._xdet = row['x']
        self._gui._sess_sel._ydet = row['y']
        x = row['x']
        xlim, ylim = self._gui._graph_det._define_lim(x)
        self._gui._graph_split = False
        self._gui._graph_det._graph._cursor_lines = []
        self._gui._graph_det._refresh(self._gui._sess_items, xlim=xlim,
                                      ylim=ylim)#, init_cursor=True)


    def _on_histogram(self, event):
        if not hasattr(self._gui, '_graph_hist'):
            from .gui_graph import GUIGraphHistogram
            GUIGraphHistogram(self._gui)#, self._col_values)
        else:
            #self._gui._graph_hist._graph._fig.clear()
            self._gui._graph_hist._fig.clear()
        self._gui._graph_hist._refresh(self._gui._sess_items)


    def _on_label_right_click(self, event):
        row, col = event.GetRow(), event.GetCol()
        if row == -1:
            self._gui._col_sel = col
            self._gui._col_tab = self._tab
            self._gui._col_values = [float(self._tab.GetCellValue(i, col)) \
                                     for i in range(self._tab.GetNumberRows())]
            self.PopupMenu(GUITablePopup(self._gui, self, event, 'Histogram',
                                         'histogram'), event.GetPosition())
        if col == -1:
            self.PopupMenu(GUITablePopup(self._gui, self, event, 'Remove',
                                         'remove'), event.GetPosition())

    def _on_remove(self, event):
        row = self._gui._tab_popup._event.GetRow()
        self._remove_data(row)
        self._gui._sess_sel.cb._spec_update()
        #self._tab.DeleteRows(pos=len(self._data.t), numRows=1)
        self._fill()
        self._gui._refresh(init_cursor=True)


    def _on_view(self, event, from_scratch=True):
        self._data = getattr(self._gui._sess_sel, self._attr)
        if 'z' in self._data.t.colnames: self._data.t.sort('z')
        if 'x' in self._data.t.colnames: self._data.t.sort('x')
        try:
            self._tab.DeleteCols(pos=0, numCols=self._tab.GetNumberCols())
            #self._tab.DeleteRows(pos=0, numRows=self._tab.GetNumberRows())
            #logging.info("I'm updating table...")
        except:
            logging.info("I'm loading table...")
        self._init(from_scratch)
        self._fill()
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        self._panel.SetSizer(self._box)
        self._tab.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self._on_detail)
        self._tab.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK,
                       self._on_label_right_click)
        self.Bind(wx.EVT_CLOSE, self._on_close)
        self.Centre()
        self.SetPosition((wx.DisplaySize()[0]*0.02, wx.DisplaySize()[1]*0.23))
        self.Show()


    def _remove_data(self, row):

        """
        if self._attr == 'systs':
            id = self._data.t['id'][row]
            for i, m in enumerate(self._data._mods_t):
                if id in m['id']:
                    if len(m['id']) == 1:
                        self._data._mods_t.remove_row(i)
                    else:
                        self._data._mods_t['id'][i].remove(id)
        """
        self._data.t.remove_row(row)
        if self._attr == 'systs':
            self._gui._sess_sel.cb._mods_update_old()


class GUITableLineList(GUITable):
    """ Class for the GUI line list """

    def __init__(self,
                 gui,
                 title="Line table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.2):

        super(GUITableLineList, self).__init__(gui, 'lines', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_lines = self

#    def _on_show(self, event):
#        pass

class GUITableModelList(GUITable):
    """ Class for the GUI model list """

    def __init__(self,
                 gui,
                 title="Model table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.2):

        super(GUITableModelList, self).__init__(gui, 'mods', title,
                                                size_x, size_y)

        self._gui = gui
        self._gui._tab_mods = self

class GUITablePopup(wx.Menu):
    """ Class for the GUI table popup menu """

    def __init__(self, gui, parent, event, title, attr):
        super(GUITablePopup, self).__init__()
        self._parent = parent
        self._event = event
        self._gui = gui
        self._gui._tab_popup = self
        self._parent = parent
        self._event = event
        self._title = np.array(title, ndmin=1)
        self._attr = np.array(attr, ndmin=1)

        #show = wx.MenuItem(self, wx.NewId(), 'Show')
        #self.Bind(wx.EVT_MENU, self._parent._on_show, show)
        #self.Append(show)
        """
        if isinstance(self._parent, GUITableSystList):
            fit = wx.MenuItem(self, wx.NewId(), 'Fit...')
            improve = wx.MenuItem(self, wx.NewId(), 'Improve...')
            self.Bind(wx.EVT_MENU, self._parent._on_fit, fit)
            self.Bind(wx.EVT_MENU, self._parent._on_improve, improve)
            self.Append(fit)
            self.Append(improve)
        remove = wx.MenuItem(self, wx.NewId(), 'Remove')
        self.Bind(wx.EVT_MENU, self._parent._on_remove, remove)
        self.Append(remove)
        """
        self._items = []
        for t, a in zip(self._title, self._attr):
            item = wx.MenuItem(self, wx.NewId(), t)
            self.Bind(wx.EVT_MENU, getattr(self._parent, '_on_'+a), item)
            self.Append(item)
            self._items.append(item)

class GUITableSpectrum(GUITable):
    """ Class for the GUI spectrum table """

    def __init__(self,
                 gui,
                 title="Spectrum table",
                 size_x=wx.DisplaySize()[0]*0.5,
                 size_y=wx.DisplaySize()[1]*0.2):

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
                 size_y=wx.DisplaySize()[1]*0.2):

        super(GUITableSystList, self).__init__(gui, 'systs', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_systs = self

    def _extract_mod(self, tab, row):
        systs = self._gui._sess_sel.systs
        labels = np.array([tab.GetColLabelValue(i).split('\n')[0] \
                  for i in range(tab.GetNumberCols())])
        id = int(tab.GetCellValue(row, np.where(labels == 'id')[0][0]))
        pars = ['z', 'logN', 'b']
        vals = [float(tab.GetCellValue(row, np.where(labels == p)[0][0])) \
                for p in pars]
        cols = [tab.GetCellTextColour(row, np.where(labels == p)[0][0]) \
                for p in pars]
        for m in systs._mods_t:
            if id in m['id']:
                for p, v, c in zip(pars, vals, cols):
                    par_name = 'lines_voigt_%i_%s' % (id, p)
                    m['mod']._pars[par_name].set(value=v)
                    m['mod']._pars[par_name].set(vary=True)
                    if c == 'grey':
                        m['mod']._pars[par_name].set(vary=False)
                return m['mod']

    def _on_cell_right_click(self, event):
        row = event.GetRow()
        col = event.GetCol()
        if col in [3, 5, 7]:
            if self._tab.GetCellTextColour(row, col) == 'black':
                title = 'Freeze parameter'
            elif self._tab.GetCellTextColour(row, col) == 'grey':
                title = 'Unfreeze parameter'
            self.PopupMenu(GUITablePopup(self._gui, self, event,
                                         title, 'freeze_par'),
                           event.GetPosition())

    def _on_detail(self, event, span=30):
        if event.GetRow() == -1: return
        if not hasattr(self._gui, '_graph_det'):
            from .gui_graph import GUIGraphDetail
            GUIGraphDetail(self._gui, init_ax=False)
        #elif len(self._gui._graph_det._graph._fig.axes) == 1:
        else:
            self._gui._graph_det._graph._fig.clear()

        # Color background of systems in the same group
        mods_sel = np.where([self._data.t['id'][event.GetRow()] in i \
                             for i in self._gui._sess_sel.systs._mods_t['id']])
        for j, r in enumerate(self._data.t):
            for i in range(len(self._data.t.colnames)):
                if r['id'] in np.array(self._gui._sess_sel.systs._mods_t['id'][mods_sel][0]):
                    self._tab.SetCellBackgroundColour(j, i, 'cyan')
                else:
                    self._tab.SetCellBackgroundColour(j, i, None)
        self._tab.ForceRefresh()


        #row = self._data.t[self._gui._tab_popup._event.GetRow()]
        row = self._data.t[event.GetRow()]
        graph = self._gui._graph_det._graph
        series = series_d[row['series']]
        rows = min(4, len(series))
        cols = len(series)//5+1
        size_x = wx.DisplaySize()[0]*0.4*cols
        size_y = min(wx.DisplaySize()[1]*0.9, wx.DisplaySize()[1]*0.3*rows)
        self._gui._graph_det.SetSize(wx.Size(size_x, size_y))

        # Redshift and wavelengths need to be initialized before the cursor
        # is created in the graph
        graph._axes = OrderedDict()
        graph._zems = OrderedDict()
        graph._xs = OrderedDict()
        graph._series = OrderedDict()
        graph._cursor_lines = []
        for i, s in enumerate(series):
            #key = s[-4:]
            key = s.split('_')[-1]
            x = (1+row['z'])*xem_d[s]
            zem = (1+row['z'])*xem_d[s]/xem_d['Ly_a']-1
            #print('out', xem_d[s], graph._zem, graph._x)
            graph._zems[key] = zem
            graph._xs[key] = x
            graph._series[key] = s
            if i == 0:
                graph._x = x
                graph._zem = zem

        #graph._axes = []
        #for i, (x, zem) in enumerate(zip(graph._xs, graph._zems)):
            xunit = self._gui._sess_sel.spec.x.unit
            self._gui._sess_sel.cb.x_convert(zem=zem)
            self._gui._sess_sel._xdet = x
            self._gui._sess_sel._ydet = 0.0
            _, ylim = self._gui._graph_det._define_lim(0)#, norm=True)

            if i == 0:
                graph._ax = graph._fig.add_subplot(rows, cols, i+1)
                title = row['series']
                if len(series) > 1:
                    graph._ax.tick_params(labelbottom=False)
            else:
                title = None
                graph._ax = graph._fig.add_subplot(rows, cols, i+1,
                                                   sharex=graph._ax,
                                                   sharey=graph._ax)
            graph._ax.tick_params(top=True, right=True, direction='in')
            graph._fig.subplots_adjust(hspace=0.)
            graph._axes[key] = graph._ax
            self._gui._graph_det._refresh(
                self._gui._sess_items, title=title, text=key,
                xlim=(-500, 500), ylim=ylim)

            self._gui._sess_sel.cb.x_convert(zem=zem, xunit=xunit)

    def _on_fit(self, event):
        row = self._gui._tab_popup._event.GetRow()
        series = self._tab.GetCellValue(row, 1)
        z = float(self._tab.GetCellValue(row, 3))
        logN = float(self._tab.GetCellValue(row, 5))
        b = float(self._tab.GetCellValue(row, 7))
        cb = self._gui._sess_sel.cb
        mod = self._extract_mod(self._tab, row)
        mod._fit()
        self._gui._sess_sel.systs._update(mod, mod_t=False)
        cb._spec_update()
        self._gui._refresh(init_cursor=True)

    def _on_freeze_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        par = self._tab.GetColLabelValue(col)[:-1]
        id = int(self._tab.GetCellValue(row, 10))
        value = float(self._tab.GetCellValue(row, col))
        if self._tab.GetCellTextColour(row, col) == 'black':
            self._tab.SetCellTextColour(row, col, 'grey')
            #self._gui._sess_sel.cb._constrain(par, id, value, False)
        elif self._tab.GetCellTextColour(row, col) == 'grey':
            self._tab.SetCellTextColour(row, col, 'black')
        self._tab.ForceRefresh()

    def _on_improve(self, event):
        row = self._gui._tab_popup._event.GetRow()
        z = float(self._tab.GetCellValue(row, 3))
        self._gui._sess_sel.add_syst_from_resids(z_start=z-1e-3, z_end=z+1e-3)
        self._gui._refresh(init_cursor=True)

    def _on_label_right_click(self, event):
        row, col = event.GetRow(), event.GetCol()
        if row == -1 and col>1:
            self._gui._col_sel = col
            self._gui._col_tab = self._tab
            self._gui._col_values = [float(self._tab.GetCellValue(i, col)) \
                                     for i in range(self._tab.GetNumberRows())]
            self.PopupMenu(GUITablePopup(self._gui, self, event, 'Histogram',
                                         'histogram'), event.GetPosition())
        if col == -1:
            self.PopupMenu(GUITablePopup(self._gui, self, event,
                                         ['Fit', 'Improve', 'Remove'],
                                         ['fit', 'improve', 'remove']),
                           event.GetPosition())

    def _on_view(self, event, **kwargs):
        super(GUITableSystList, self)._on_view(event, **kwargs)
        self._tab.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK,
                       self._on_cell_right_click)
