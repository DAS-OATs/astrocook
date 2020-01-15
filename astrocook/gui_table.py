from .functions import get_selected_cells, trans_parse
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
                    self._tab.SetCellValue(j, i, "%3.7f" % r[n])
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


    def _labels_extract(self):
        return np.array([self._tab.GetColLabelValue(i).split('\n')[0] \
                         for i in range(self._tab.GetNumberCols())])


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


    def _on_edit(self, event):
        row, col = event.GetRow(), event.GetCol()
        labels = self._labels_extract()
        self._data.t[labels[col]][row] = self._tab.GetCellValue(row, col)


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
            self._gui._col_values = [#float(self._tab.GetCellValue(i, col)) \
                                     self._data.t[self._labels_extract()[col]][i] \
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
        self._tab.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self._on_edit)
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
        if self._attr == 'systs':
            #self._gui._sess_sel.cb._mods_update_old()
            #self._gui._sess_sel.cb._mods_recreate()
            self._gui._sess_sel.cb._systs_remove([row])
            #self._gui._sess_sel.cb._systs_refit(refit_id, max_nfev_def)
            self._gui._sess_sel.cb._systs_cycle()
        else:
            self._data.t.remove_row(row)




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
            if t == 'sep':
                self.AppendSeparator()
            else:
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
        self._freezes_l = []
        self._links_l = []
        self._freezes_d = {}
        self._links_d = {}


    def _constr_copy(self):
        labels = self._labels_extract()
        for i in range(self._tab.GetNumberRows()):
            id = self._id_extract(i)
            for m in self._gui._sess_sel.systs._mods_t:
                if id in m['id']:
                    for p,v in m['mod']._pars.items():
                        #"""
                        if v.vary == False:
                            try:
                                col = np.where(labels==p.split('_')[-1])[0][0]
                                id_i = i if col==9 else np.where(id==int(p.split('_')[-2]))[0][0]
                                self._tab.SetCellTextColour(i, col, 'grey')
                            except:
                                pass
                        #"""
                        if v.expr is not None:
                            try:
                                col = np.where(np.logical_or(
                                    labels==p.split('_')[-1],
                                    labels==v.expr.split('_')[-1]))[0]
                                for c in col:
                                    self._tab.SetCellTextColour(i, c,
                                                                'forest green')
                            except:
                                pass

    def _id_extract(self, row):
        labels = self._labels_extract()
        return int(self._tab.GetCellValue(row, np.where(labels == 'id')[0][0]))


    def _key_extract(self, row, col):
        #id = self._tab.GetCellValue(row, 11).strip()
        id = self._id_extract(row)
        coln = self._tab.GetColLabelValue(col).split('\n')[0]
        parn = 'lines_%s_%s_%s' % (self._tab.GetCellValue(row, 0), id, coln)
        return int(id), parn


    def _mod_extract(self, row):
        """
        systs = self._gui._sess_sel.systs
        labels = np.array([self._tab.GetColLabelValue(i).split('\n')[0] \
                  for i in range(self._tab.GetNumberCols())])
        id = int(tab.GetCellValue(row, np.where(labels == 'id')[0][0]))
        """
        id = self._id_extract(row)
        """
        pars = ['z', 'logN', 'b']
        vals = [float(self._tab.GetCellValue(row, np.where(labels == p)[0][0])) \
                for p in pars]
        cols = [self._tab.GetCellTextColour(row, np.where(labels == p)[0][0]) \
                for p in pars]
        """
        for m in self._gui._sess_sel.systs._mods_t:
            if id in m['id']:
                return m['mod']


    def _mods_edit(self, dict):
        for k, v in dict.items():
            #print(k, dict[k])
            for m in self._gui._sess_sel.systs._mods_t:
                if v[0] in m['id']:
                    if v[1]=='expr': m['mod']._pars[k].set(expr=v[2])
                    if v[1]=='vary': m['mod']._pars[k].set(vary=v[2])

    def _on_cell_right_click(self, event):
        row = event.GetRow()
        col = event.GetCol()
        sel = get_selected_cells(self._tab)
        if len(sel) == 1:
            self._tab.SetGridCursor(row, col)
            sel = get_selected_cells(self._tab)
        self._cells_sel = []
        for s in sel:
            if s[1] in [3, 5, 7]: self._cells_sel.append(s)
        if col in [3, 5, 7]:
            title = []
            attr = []
            if len(self._cells_sel) > 1:
                if self._tab.GetCellTextColour(row, col) == 'forest green':
                    title = ['Unlink']
                else:
                    title = ['Link']
                attr = ['link_par']
            if self._tab.GetCellTextColour(row, col) == 'grey':
                title += ['Unfreeze']
            else:
                title += ['Freeze']
            attr += ['freeze_par']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
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
        #print(self._gui._sess_sel.systs._mods_t['id'])
        #print(mods_sel)
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
        try:
            series = series_d[row['series']]
        except:
            series = row['series'].split(',')
        series = trans_parse(row['series'])
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


    def _on_edit(self, event):
        row, col = event.GetRow(), event.GetCol()
        value = float(self._tab.GetCellValue(row, col))
        labels = self._labels_extract()
        self._data.t[labels[col]][row] = value
        id = self._id_extract(row)
        mod = self._mod_extract(row)
        mod._pars['lines_voigt_%i_%s' % (id, labels[col])].set(value=value)


    def _on_ccf(self, event):
        row = self._gui._tab_popup._event.GetRow()
        cb = self._gui._sess_sel.cb
        mod = self._mod_extract(row)
        cb._mod_ccf(mod)


    def _on_ccf_max(self, event):
        row = self._gui._tab_popup._event.GetRow()
        cb = self._gui._sess_sel.cb
        mod = self._mod_extract(row)
        cb._mod_ccf_max(mod)


    def _on_fit(self, event):
        row = self._gui._tab_popup._event.GetRow()
        series = self._tab.GetCellValue(row, 1)
        z = float(self._tab.GetCellValue(row, 3))
        logN = float(self._tab.GetCellValue(row, 5))
        b = float(self._tab.GetCellValue(row, 7))
        cb = self._gui._sess_sel.cb
        mod = self._mod_extract(row)
        #mod._fit()
        cb._syst_fit(mod, max_nfev_def)
        #self._gui._sess_sel.systs._update(mod, mod_t=False)
        #cb._systs_refit(max_nfev=max_nfev_def)
        #cb._systs_cycle()
        cb._spec_update()
        self._gui._refresh(init_cursor=True)



    def _on_freeze_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        for (r, c) in self._cells_sel:
            id, parn = self._key_extract(r, c)
            if self._tab.GetCellTextColour(row, col) == 'grey':
                self._tab.SetCellTextColour(r, c, 'black')
                try:
                    self._freezes_l.remove((r,c))
                    del self._freezes_d[parn]
                except:
                    pass
            else:
                #self._tab.SetCellTextColour(r, c, 'grey')
                self._freezes_l.append((r,c))
                self._freezes_d[parn] = (id, 'vary', False)
        self._tab.ForceRefresh()
        self._mods_edit(self._freezes_d)
        self._constr_copy()


    def _on_freeze_par_all(self, event=None, col=None):
        if event is not None:
            col = self._gui._tab_popup._event.GetCol()
        for i in range(self._tab.GetNumberRows()):
            id, parn = self._key_extract(i, col)
            #self._tab.SetCellTextColour(i, col, 'grey')
            self._freezes_l.append((i, col))
            self._freezes_d[parn] = (id, 'vary', False)
        self._tab.ForceRefresh()
        self._mods_edit(self._freezes_d)
        self._constr_copy()


    def _on_improve(self, event):
        row = self._gui._tab_popup._event.GetRow()
        z = float(self._tab.GetCellValue(row, 3))
        self._gui._sess_sel.cb.systs_new_from_resids(z_start=z-1e-3, z_end=z+1e-3)
        self._gui._refresh(init_cursor=True)


    def _on_label_right_click(self, event):
        row, col = event.GetRow(), event.GetCol()
        if row == -1 and col>1:
            self._gui._col_sel = col
            self._gui._col_tab = self._tab
            self._gui._col_values = [#float(self._tab.GetCellValue(i, col)) \
                                     self._data.t[self._labels_extract()[col]][i] \
                                     for i in range(self._tab.GetNumberRows())]
            title = ['Histogram']
            attr = ['histogram']
            if col in [3,5,7]:
                title += ['Freeze all']
                attr += ['freeze_par_all']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                           event.GetPosition())
        if col == -1:
            self.PopupMenu(GUITablePopup(
                self._gui, self, event,
                ['Fit', 'Improve', 'Remove', 'sep', 'CCF', 'Maximize CCF'],
                ['fit', 'improve', 'remove', None, 'ccf', 'ccf_max']),
                event.GetPosition())


    def _on_link_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        for (r, c) in self._cells_sel:
            """
            key = 'lines_%s_%s_%s' % (self._tab.GetCellValue(r, 0),
                                      self._tab.GetCellValue(r, 11).strip(),
                                      self._tab.GetColLabelValue(c).split('\n')[0])
            val = 'lines_%s_%s_%s' % (self._tab.GetCellValue(row, 0),
                                      self._tab.GetCellValue(row, 11).strip(),
                                      self._tab.GetColLabelValue(col).split('\n')[0])
            """
            id_1, parn_1 = self._key_extract(r, c)
            id_2, parn_2 = self._key_extract(row, col)
            if id_1<id_2: id, parn, val = id_2, parn_2, parn_1
            if id_1>=id_2: id, parn, val = id_1, parn_1, parn_2
            if self._tab.GetCellTextColour(row, col) == 'forest green':
                self._tab.SetCellTextColour(r, c, 'black')
                if parn != val:
                    try:
                        self._links_l.remove((r,c))
                        del self._links_d[parn]
                    except:
                        pass
            else:
                self._tab.SetCellTextColour(r, c, 'forest green')
                if parn != val:
                    self._links_l.append((r,c))
                    self._links_d[parn] = (id, 'expr', val)
        self._tab.ForceRefresh()
        self._constr_copy()
        self._mods_edit(self._links_d)


    def _on_view(self, event, **kwargs):
        super(GUITableSystList, self)._on_view(event, **kwargs)
        self._constr_copy()
        #self._tab.Bind(wx.grid.EVT_GRID_CELL_LEFT_DCLICK,
        #               self._on_cell_left_dclick)
        self._tab.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK,
                       self._on_cell_right_click)
