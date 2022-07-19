from .functions import get_selected_cells, to_x, to_z, trans_parse
from .gui_dialog import *
from .vars import *
from collections import OrderedDict
import cProfile
import logging
import matplotlib.pyplot as plt
import numpy as np
import pprint
import pstats
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix
import wx.lib.colourdb as cdb

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
        self._tab_id = self._gui._menu_tab_id
        self._gui._tab = self
        self._menu = self._gui._panel_sess._menu._view._menu
        #self._open = {'spec': False, 'lines': False, 'systs': False}
        super(GUITable, self).__init__(parent=None, title=self._title,
                                       size=(self._size_x, self._size_y))
        self._shown = False


    def _data_edit(self, row, label, value, attr=None):
        if attr is None: attr = self._attr
        tab = getattr(self._gui, '_tab_'+attr)
        tab._data.t[label][row] = value


    def _data_init(self, from_scratch=True, autosort=False, attr=None):
        if attr is None: attr = self._attr
        tab = getattr(self._gui, '_tab_'+attr)
        sess = self._gui._sess_sel
        tab._data = getattr(sess, attr)
        if autosort:
            if 'z' in tab._data.t.colnames: tab._data.t.sort(['z','id'])
            if 'x' in tab._data.t.colnames: tab._data.t.sort('x')
        try:
            tab._tab.DeleteCols(pos=0, numCols=tab._tab.GetNumberCols())
        except:
            logging.info("I'm loading table...")
        self._init(from_scratch, attr)
        self._fill(attr)


    def _data_remove(self, row, attr=None, log=True):
        sess = self._gui._sess_sel
        if attr is None: attr = self._attr

        if self._attr == 'systs':
            sess.cb._systs_remove([row])
            sess.cb._mods_recreate()
            sess.cb._spec_update()
        else:
            tab = getattr(self._gui, '_tab_'+attr)
            tab._data.t.remove_row(row)


    def _data_sort(self, label, reverse=False, attr=None, log=True):
        if attr is None: attr = self._attr
        tab = getattr(self._gui, '_tab_'+attr)

        if reverse and self._attr=="systs":
            tab._data.t['id'] = -1*tab._data.t['id']
        if self._attr=="systs":
            tab._data.t.sort([label, 'id'], reverse=reverse)
        else:
            tab._data.t.sort(label, reverse=reverse)
        if reverse and self._attr=="systs":
            tab._data.t['id'] = -1*tab._data.t['id']


    def _fill(self, attr=None):
        if attr is None: attr = self._attr

        tab = getattr(self._gui, '_tab_'+attr)
        for j, r in enumerate(tab._data.t):
            for i, n in enumerate(tab._data.t.colnames):

                if j == 0:
                    tab._tab.SetColSize(i, 150)
                    tab._tab.SetColLabelValue(i, "%s\n%s" \
                                              % (n, str(tab._data.t[n].unit)))
                if type(r[n]) == np.int64:
                    tab._tab.SetCellValue(j, i, "%4i" % r[n])
                elif type(r[n]) == str or type(r[n]) == np.str_:
                    tab._tab.SetCellValue(j, i, r[n])
                elif type(r[n]) == OrderedDict:
                    tab._tab.SetCellValue(j, i, pprint.pformat(r[n]))
                elif type(r[n]) == dict:
                    tab._tab.SetCellValue(j, i, pprint.pformat(r[n]))
                else:
                    if n in ['logN', 'dlogN', 'b', 'db', 'btur', 'dbtur', 'resol', 'chi2r', \
                             'snr']:
                        format = '%3.3'
                    else:
                        format = '%3.7'
                    if np.abs(r[n])<1e-7 and r[n]!=0:
                        #tab._tab.SetCellValue(j, i, "%3.7e" % r[n])
                        format += 'e'
                    else:
                        #tab._tab.SetCellValue(j, i, "%3.7f" % r[n])
                        format += 'f'
                    tab._tab.SetCellValue(j, i, format % r[n])
        tab._tab.AutoSizeColumns(True)


    def _init(self, from_scratch=True, attr=None):
        if attr is None: attr = self._attr

        tab = getattr(self._gui, '_tab_'+attr)

        if not from_scratch:
            try:
                if tab._tab.GetNumberRows() != 0:
                    tab._tab.DeleteRows(0, tab._tab.GetNumberRows())
            except:
                from_scratch = True
        if from_scratch:
            super(GUITable, tab).__init__(parent=None, title=tab._title,
                                           size=(tab._size_x, tab._size_y))

            tab._panel = wx.Panel(tab)
            tab._tab = gridlib.Grid(tab._panel)
            tab._tab.CreateGrid(0, 0)
            tab.SetPosition((0, wx.DisplaySize()[1]*0.5))

        coln = len(tab._data.t.colnames)
        rown = len(tab._data.t)-tab._tab.GetNumberRows()
        tab._tab.AppendCols(coln)
        tab._tab.AppendRows(rown)


    def _labels_extract(self):
        if not hasattr(self, '_tab'):
            self._init(False)
        return np.array([self._tab.GetColLabelValue(i).split('\n')[0] \
                         for i in range(self._tab.GetNumberCols())])


    def _on_close(self, event):
        self._shown = False
        self.Destroy()
        #self.Close()


    def _on_detail(self, event):
        if self._attr != 'systs': return
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
        #self._gui._graph_det._graph._cursor_lines = []
        self._gui._graph_det._refresh(self._gui._sess_items, xlim=xlim,
                                      ylim=ylim)#, init_cursor=True)


    def _on_edit(self, event):
        row, col = event.GetRow(), event.GetCol()
        labels = self._labels_extract()
        label = labels[col]
        value = self._tab.GetCellValue(row, col)
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab", "_data_init",
                                            {"attr": self._attr})
        sess.json += self._gui._json_update("_tab", "_data_edit",
                                            {"row": row, "label": label,
                                             "value": value,
                                             "attr": self._attr})
        """
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab', '_data_edit',
                             {'row': row, 'label': label, 'value': value,
                              'attr': self._attr})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_edit(row, label, value, self._attr)


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
            self._gui._col_unit = self._data.t[self._labels_extract()[col]].unit
            title = ['Sort ascending', 'Sort descending', 'sep', 'Histogram']
            attr = ['sort', 'sort_reverse', None, 'histogram']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title,
                                         attr), event.GetPosition())
        if col == -1:
            self.PopupMenu(GUITablePopup(self._gui, self, event, 'Remove',
                                         'remove'), event.GetPosition())

    def _on_remove(self, event):
        row = self._gui._tab_popup._event.GetRow()

        sess = self._gui._sess_sel
        if self._attr == 'systs':
            """
            sess.json += self._gui._json_update("_tab", "_data_init",
                                                {"attr": self._attr})
            sess.json += self._gui._json_update("cb", "_systs_remove",
                                                {"rem": [row]})
            sess.json += self._gui._json_update("cb", "_mods_recreate", {})
            """
            sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
            sess.log.append_full('cb', '_systs_remove', {'rem': [row]})
            sess.log.append_full('cb', '_mods_recreate', {})
        else:
            """
            sess.json += self._gui._json_update("_tab", "_data_init",
                                                {"attr": self._attr})
            sess.json += self._gui._json_update("_tab", "_data_remove",
                                                {"row": row, "attr": self._attr})
            """
            sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
            sess.log.append_full('_tab', '_data_remove',
                                 {'row': row, 'attr': self._attr})
        """
        sess.json += self._gui._json_update("cb", "_spec_update", {})
        """
        sess.log.append_full('cb', '_spec_update', {})

        self._data_remove(row, self._attr)
        #sess.cb._spec_update()
        #self._tab.DeleteRows(pos=len(self._data.t), numRows=1)
        #self._fill()
        self._gui._refresh(init_cursor=True)


    def _on_sort(self, event):
        labels = self._labels_extract()
        #self._data.t.sort([labels[self._gui._col_sel], 'id'])

        label = labels[self._gui._col_sel]
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab", "_data_init",
                                            {"attr": self._attr})
        sess.json += self._gui._json_update("_tab", "_data_sort",
                                            {"label": label,
                                             "attr": self._attr})
        """
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab', '_data_sort',
                             {'label': label, 'attr': self._attr})
        self._data_sort(label=label, attr=self._attr)
        self._gui._refresh(autosort=False)


    def _on_sort_reverse(self, event):
        labels = self._labels_extract()
        label = labels[self._gui._col_sel]
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab", "_data_init",
                                            {"attr": self._attr})
        sess.json += self._gui._json_update("_tab", "_data_sort",
                                            {"label": label,
                                             "attr": self._attr,
                                             "reverse": True})
        """
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab', '_data_sort',
            {'label': label, 'attr': self._attr, 'reverse': True})
        self._data_sort(label=label, attr=self._attr, reverse=True)
        self._gui._refresh(autosort=False)


    def _on_view(self, event=None, from_scratch=True, autosort=False):
        #sess = self._gui._sess_sel
        #print('on_view')
        #sess.json += self._gui._json_update("_tab", "_data_init",
        #                                    {"attr": self._attr})
        #sess.log.append_full('_tab', '_data_init', {'attr': self._attr})
        #if hasattr(self, '_dlg_mini_log') and self._dlg_mini_log._shown:
        #    self._gui._dlg_mini_log._refresh()
        self._view(event, from_scratch, autosort)

    def _view(self, event=None, from_scratch=True, autosort=False):
        self._data_init(from_scratch, autosort)
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
        self._shown = True
        #print(self, self._shown)


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


    def _on_close(self, event, **kwargs):
        super(GUITableLineList, self)._on_close(event, **kwargs)
        #del self._gui._tab_lines._data
        self._menu.FindItemById(self._tab_id[1]).Check(False)


    def _on_view(self, event, **kwargs):
        super(GUITableLineList, self)._on_view(event, **kwargs)


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
        self._gui = gui
        self._gui._tab_popup = self
        self._parent = parent
        self._event = event
        self._title = np.array(title, ndmin=1)
        self._attr = np.array(attr, ndmin=1)

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


    def _on_close(self, event, **kwargs):
        super(GUITableSpectrum, self)._on_close(event, **kwargs)
        #del self._gui._tab_spec._data
        self._menu.FindItemById(self._tab_id[0]).Check(False)


    def _on_view(self, event, **kwargs):
        super(GUITableSpectrum, self)._on_view(event, **kwargs)


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
        self._freezes_d = {}
        self._links_d = {}
        self._colours = ['cadet blue', 'forest green', 'dark orchid', 'purple', 'maroon']#, 'orchid']
        self._colourc = 0
        self._links_c = {}
        self._cells_sel = []

    def _data_cell_right_click(self, row, col):
        sel = get_selected_cells(self._tab)
        if len(sel) == 1:
            self._tab.SetGridCursor(row, col)
            sel = get_selected_cells(self._tab)
        #self._cells_sel = []
        self._data_cells_desel()
        for s in sel:
            if s[1] in [3, 5, 7, 9]: #self._cells_sel.append(s)
                row = s[0]
                col = s[1]
                self._data_cells_sel(row, col)


    def _data_cells_sel(self, row, col, log=True):
        """
        if json:
            sess = self._gui._sess_sel
            sess.json += self._gui._json_update("_tab_systs", "_data_cells_sel",
                                                {"row": row, "col": col,
                                                 "json": False})
        """
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
            sess.log.append_full('_tab_systs', '_data_cells_sel',
                                 {'row': row, 'col': col, 'log': False})
            if hasattr(self._gui, '_dlg_mini_log') \
                and self._gui._dlg_mini_log._shown:
                self._gui._dlg_mini_log._refresh()
        self._cells_sel.append((row,col))


    def _data_cells_desel(self, log=True):
        """
        if json:
            sess = self._gui._sess_sel
            sess.json += self._gui._json_update("_tab_systs", "_data_cells_desel",
                                                {"json": False})
        """
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
            sess.log.append_full('_tab_systs', '_data_cells_desel',
                                 {'log': False})
            if hasattr(self._gui, '_dlg_mini_log') \
                and self._gui._dlg_mini_log._shown:
                self._gui._dlg_mini_log._refresh()
        self._cells_sel = []


    def _data_detail(self, row_z, row_series, row_id, span=30):
        if not hasattr(self._gui, '_graph_det'):
            from .gui_graph import GUIGraphDetail
            GUIGraphDetail(self._gui, init_ax=False)
        else:
            self._gui._graph_det._graph._fig.clear()

        z = row_z #row['z']
        series = trans_parse(row_series) #row['series'])
        self._gui._graph_main._z_sel = z
        self._gui._graph_main._series_sel = series
        self._gui._graph_main._refresh(self._gui._sess_sel)
        self._gui._graph_det._update(series, z, hwin_def)
        if not hasattr(self._gui, '_dlg_mini_systems') \
            or self._gui._dlg_mini_systems == None:
            #GUIDialogMiniSystems(self._gui, "System controls", series=row['series'], z=row['z'])
            GUIDialogMiniSystems(self._gui, "System controls", series=row_series, z=row_z)
        else:
            #self._gui._dlg_mini_systems._refresh(row['series'], row['z'])
            self._gui._dlg_mini_systems._refresh(row_series, row_z)

        # Color background of systems in the same group
        #mods_sel = np.where([self._data.t['id'][event.GetRow()] in i \
        mods_sel = np.where([row_id in i \
                             for i in self._gui._sess_sel.systs._mods_t['id']])
        for j, r in enumerate(self._data.t):
            for i in range(len(self._data.t.colnames)):
                if r['id'] in np.array(self._gui._sess_sel.systs._mods_t['id'][mods_sel][0]):
                    self._tab.SetCellBackgroundColour(j, i, 'cyan')
                else:
                    self._tab.SetCellBackgroundColour(j, i, None)
        self._tab.ForceRefresh()

    def _data_edit(self, row, label, value, update_mod=True):
        self._data.t[label][row] = value
        if update_mod:
            id = self._id_extract(row)
            mod = self._mod_extract(row)
            if label == 'series':
                z0 = float(self._tab.GetCellValue(row, 3))
                z = mod._pars['lines_voigt_%i_z' % id].value
                z0_new = to_z(to_x(z0, trans_parse(mod._series)[0]),
                             trans_parse(value)[0])
                z_new = to_z(to_x(z, trans_parse(mod._series)[0]),
                             trans_parse(value)[0])
                vary = mod._pars['lines_voigt_%i_z' % id].vary
                expr = mod._pars['lines_voigt_%i_z' % id].expr
                mod._pars['lines_voigt_%i_z' % id].set(
                    value=z_new, vary=vary, expr=expr)
                self._tab.SetCellValue(row, 2, "%.7f" % z0_new)
                self._tab.SetCellValue(row, 3, "%.7f" % z_new)
                self._data.t['z0'][row] = z0_new
                self._data.t['z'][row] = z_new
                mod._series = value
                self._gui._sess_sel.cb._mods_recreate2()
            else:
                try:
                    vary = mod._pars['lines_voigt_%i_%s' % (id, label)].vary
                    expr = mod._pars['lines_voigt_%i_%s' % (id, label)].expr
                    mod._pars['lines_voigt_%i_%s' % (id, label)].set(
                        value=value, vary=vary, expr=expr)
                except:
                    vary = mod._pars['psf_gauss_%i_%s' % (id, label)].vary
                    expr = mod._pars['psf_gauss_%i_%s' % (id, label)].expr
                    mod._pars['psf_gauss_%i_%s' % (id, label)].set(
                        value=value, vary=vary, expr=expr)


    def _data_fit(self, row):
        cb = self._gui._sess_sel.cb
        mod = self._mod_extract(row)
        cb._syst_fit(mod, max_nfev_def)
        cb._spec_update()


    def _data_freeze_par(self, row, col):
        for (r, c) in self._cells_sel:
            id, parn = self._key_extract(r, c)
            if self._tab.GetCellTextColour(row, col) != 'black':
                self._freezes_d[parn] = (id, 'vary', True)
            else:
                self._freezes_d[parn] = (id, 'vary', False)
        self._tab.ForceRefresh()
        systs = self._gui._sess_sel.systs

        for v in self._freezes_d:
            if v in self._links_d and self._links_d[v][2] != '' and self._freezes_d[v][2] == True:
                self._freezes_d[v] = (self._freezes_d[v][0],
                                      self._freezes_d[v][1], False)
        systs._constrain(self._freezes_d)
        self._text_colours()


    def _data_freeze_par_all(self, col):
        for i in range(self._tab.GetNumberRows()):
            id, parn = self._key_extract(i, col)
            self._freezes_d[parn] = (id, 'vary', False)
        self._tab.ForceRefresh()
        self._gui._sess_sel.systs._constrain(self._freezes_d)
        self._text_colours()

    def _data_init(self, from_scratch=True, autosort=False, attr=None):
        super(GUITableSystList, self)._data_init(from_scratch, autosort, attr)
        labels = self._labels_extract()
        self._ids = np.array([int(float(self._tab.GetCellValue(
                              i, np.where(labels == 'id')[0][0]))) \
                              for i in range(self._tab.GetNumberRows())])

    def _data_link_bt(self, row, col):
        self._cells_sel = sorted(self._cells_sel, key=lambda tup: tup[0])
        ref = np.argmin([int(self._key_extract(r,c)[0]) for r,c in self._cells_sel])
        others = [self._cells_sel[c] for c in np.setdiff1d(range(len(self._cells_sel)), [ref])]
        id_r, val = self._key_extract(self._cells_sel[ref][0], self._cells_sel[ref][1])
        labels = self._labels_extract()

        trans_r = self._trans_extract(self._cells_sel[ref][0])
        if not self._trans_check(trans_r): return 0
        mass_r = mass_d[trans_r[0]]

        for (r, c) in others:
            id, parn = self._key_extract(r, c)
            if self._tab.GetCellTextColour(r, c) in self._colours:#== 'forest green':
                self._links_d[parn] = (id, 'expr', '')
                self._links_d[val] = (id_r, 'expr', '')
            else:
                if parn != val:
                    trans = self._trans_extract(r)
                    if not self._trans_check(trans): return 0
                    mass = mass_d[trans[0]]

                    v = "%.3f" % (float(self._tab.GetCellValue(self._cells_sel[ref][0],
                                                         self._cells_sel[ref][1])) \
                                        * np.sqrt(mass_r/mass))

                    self._tab.SetCellValue(r, c, v)
                    self._links_d[parn] = (id, 'expr',
                        val+'*%.14f' % (np.sqrt(mass_r/mass)))
                    self._data_edit(r, labels[c], v, update_mod=False)
        self._tab.ForceRefresh()
        systs = self._gui._sess_sel.systs
        systs._constrain(self._links_d)
        #print(self._links_d)
        self._gui._sess_sel.cb._mods_recreate2(only_constr=True)
        self._text_colours()


    def _data_link_par(self, row, col):
        self._cells_sel = sorted(self._cells_sel, key=lambda tup: tup[0])
        ref = np.argmin([int(self._key_extract(r,c)[0]) for r,c in self._cells_sel])
        others = [self._cells_sel[c] for c in np.setdiff1d(range(len(self._cells_sel)), [ref])]
        id_r, val = self._key_extract(self._cells_sel[ref][0], self._cells_sel[ref][1])
        labels = self._labels_extract()
        for (r, c) in others:
            id, parn = self._key_extract(r, c)
            if self._tab.GetCellTextColour(r, c) in self._colours:#== 'forest green':
                self._links_d[parn] = (id, 'expr', '')
                self._links_d[val] = (id_r, 'expr', '')
            else:
                if parn != val:
                    v = self._tab.GetCellValue(self._cells_sel[ref][0],
                                               self._cells_sel[ref][1])
                    self._tab.SetCellValue(r, c, v)
                    self._links_d[parn] = (id, 'expr', val)
                    self._data_edit(r, labels[c], v, update_mod=False)
        self._tab.ForceRefresh()
        systs = self._gui._sess_sel.systs
        #print(self._freezes_d)
        #print(self._links_d)
        systs._constrain(self._links_d)
        self._gui._sess_sel.cb._mods_recreate2(only_constr=True)
        self._text_colours()


    def _data_top_label_right_click(self, col):
        self._gui._col_sel = col
        self._gui._col_tab = self._tab
        self._gui._col_values = [#float(self._tab.GetCellValue(i, col)) \
                                 self._data.t[self._labels_extract()[col]][i] \
                                 for i in range(self._tab.GetNumberRows())]


    def _id_extract(self, row):
        labels = self._labels_extract()
        return int(float(self._tab.GetCellValue(
                   row, np.where(labels == 'id')[0][0])))


    def _key_extract(self, row, col):
        id = self._id_extract(row)
        coln = self._tab.GetColLabelValue(col).split('\n')[0]
        parn = 'lines_%s_%s_%s' % (self._tab.GetCellValue(row, 0), id, coln)
        return int(id), parn


    def _mod_extract(self, row):
        id = self._id_extract(row)
        for m in self._gui._sess_sel.systs._mods_t:
            if id in m['id']:
                return m['mod']

    def _trans_check(self, trans):
        check = len(np.unique([t.split('_')[0] for t in trans]))==1
        if not check:
            logging.error("I cannot link temperatures if the series contains "
                          "different ions")
        return check

    def _trans_extract(self, row):
        return trans_parse(self._tab.GetCellValue(row, 1))


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


    def _on_cell_right_click(self, event):
        row = event.GetRow()
        col = event.GetCol()
        self._data_cell_right_click(row, col)
        if col in [3, 5, 7, 9]:
            title = []
            attr = []
            if len(self._cells_sel) > 1:
                if self._tab.GetCellTextColour(row, col) in self._colours:
                    title = ['Unlink']
                else:
                    title = ['Link']
                attr = ['link_par']
                if col == 7:
                    if self._tab.GetCellTextColour(row, col) in self._colours:
                        title = ['Unlink']
                    else:
                        title += ['Link temperatures']
                    attr += ['link_bt']
            if self._tab.GetCellTextColour(row, col) == 'grey':
                title += ['Unfreeze']
            else:
                title += ['Freeze']
            attr += ['freeze_par']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                           event.GetPosition())


    def _on_close(self, event, **kwargs):
        super(GUITableSystList, self)._on_close(event, **kwargs)
        #print(self)
        #self.Destroy()
        #self.Close()
        #del self._gui._tab_systs._data
        self._menu.FindItemById(self._tab_id[2]).Check(False)
        #self._open['systs'] = False


    def _on_detail(self, event, span=30, log=True):
        if event.GetRow() == -1: return
        row = self._data.t[event.GetRow()]
        row_z = row['z']
        row_series = row['series']
        row_id = int(row['id'])
        sess = self._gui._sess_sel
        sess.log.append_full('_tab_systs', '_data_detail',
                             {'row_z': row_z, 'row_series': row_series,
                              'row_id': row_id, 'span': span})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_detail(row_z, row_series, row_id, span)


    def _on_edit(self, event):
        row, col = event.GetRow(), event.GetCol()
        labels = self._labels_extract()
        label = labels[col]
        if col>=1:
            value = str(self._tab.GetCellValue(row, col))
        else:
            value = float(self._tab.GetCellValue(row, col))
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab_systs", "_data_edit",
                                            {"row": row, "label": label,
                                             "value": value})
        """
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab_systs', '_data_edit',
                             {'row': row, 'label': label, 'value': value})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_edit(row, label, value)


    def _on_fit(self, event):
        row = self._gui._tab_popup._event.GetRow()
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab_systs", "_data_fit",
                                            {"row": row})
        """
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab_systs', '_data_fit', {'row': row})
        self._data_init(from_scratch=False, attr='systs')
        self._data_fit(row)
        self._gui._refresh(init_cursor=True)


    def _on_fit_dialog(self, event):
        row = self._data.t[self._gui._tab_popup._event.GetRow()]
        #print(row._index)
        #params = [{'series': row['series'], 'z': "%3.7f" % float(row['z']),
        #           'logN': row['logN'], 'b': row['b'], 'refit_n': 0}]
        dlg = GUIDialogMethod(self._gui, 'Fit system', 'syst_fit',
                              params_last=[{'num': row._index+1}])

        self._gui._refresh(init_cursor=True)


    def _on_freeze_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab_systs",
                                            "_data_freeze_par",
                                            {"row": row, "col": col})
        """
        sess.log.append_full('_tab_systs', '_data_freeze_par',
                             {'row': row, 'col': col})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_freeze_par(row, col)


    def _on_freeze_par_all(self, event=None, col=None):
        if event is not None:
            col = self._gui._tab_popup._event.GetCol()
            sess = self._gui._sess_sel
            """
            sess.json += self._gui._json_update("_tab_systs",
                                                "_data_freeze_par_all",
                                                {"col": col})
            """
            sess.log.append_full('_tab_systs', '_data_freeze_par_all',
                                 {'col': col})
            if hasattr(self._gui, '_dlg_mini_log') \
                and self._gui._dlg_mini_log._shown:
                self._gui._dlg_mini_log._refresh()
            self._data_freeze_par_all(col)


    def _on_improve(self, event):
        #row = self._gui._tab_popup._event.GetRow()
        #z = float(self._tab.GetCellValue(row, 3))
        """
        sess = self._gui._sess_sel
        sess.json += self._gui._json_update("cb", "systs_improve", {})
        self._gui._sess_sel.cb.systs_improve()
        """
        dlg = GUIDialogMethod(self._gui, 'Improve systems', 'systs_improve')
        self._gui._refresh(init_cursor=True)


    def _on_label_right_click(self, event):
        row, col = event.GetRow(), event.GetCol()
        if row == -1 and col>1:
            """
            self._gui._col_sel = col
            self._gui._col_tab = self._tab
            self._gui._col_values = [#float(self._tab.GetCellValue(i, col)) \
                                     self._data.t[self._labels_extract()[col]][i] \
                                     for i in range(self._tab.GetNumberRows())]
            """
            self._data_top_label_right_click(col)
            title = ['Sort ascending', 'Sort descending', 'sep', 'Histogram']
            attr = ['sort', 'sort_reverse', None, 'histogram']
            if col in [3,5,7]:
                title += ['sep', 'Freeze all']
                attr += [None, 'freeze_par_all']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                           event.GetPosition())
        if col == -1:
            if self._gui._sess_sel.systs._compressed:
                title = ['Fit', 'Fit (open dialog)', 'Remove', 'Merge', 'sep',
                         'CCF', 'Maximize CCF', 'sep', 'Improve all']
                attr = ['fit', 'fit_dialog', 'remove', 'merge', None, 'ccf',
                        'ccf_max', None, 'improve']
            else:
                title = ['Fit', 'Fit (open dialog)', 'Remove', 'sep', 'CCF',
                         'Maximize CCF', 'sep', 'Improve all']
                attr = ['fit', 'fit_dialog', 'remove', None, 'ccf', 'ccf_max',
                        None, 'improve']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                event.GetPosition())


    def _on_link_bt(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab_systs", "_data_link_bt",
                                            {"row": row, "col": col})
        """
        sess.log.append_full('_tab_systs', '_data_link_bt',
                             {'row': row, 'col': col})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_link_bt(row, col)


    def _on_link_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        sess = self._gui._sess_sel
        """
        sess.json += self._gui._json_update("_tab_systs", "_data_link_par",
                                            {"row": row, "col": col})
        """
        sess.log.append_full('_tab_systs', '_data_link_par',
                             {'row': row, 'col': col})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_link_par(row, col)


    def _on_merge(self, event):
        row = self._data.t[self._gui._tab_popup._event.GetRow()]
        dlg = GUIDialogMethod(self._gui, 'Merge systems', 'systs_merge',
                              params_last=[{'num1': row._index+1}])

        self._gui._refresh(init_cursor=True)

    def _on_view(self, event, **kwargs):
        profile = cProfile.Profile()
        profile.enable()

        super(GUITableSystList, self)._on_view(event, **kwargs)
        #self._open['systs'] = True
        #sess = self._gui._sess_sel
        #sess.json += self._gui._json_update("_tab", "_data_init",
        #                                    {"attr": "systs"})
        self._text_colours()
        self._tab.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK,
                       self._on_cell_right_click)
        for k, v in self._data._constr.items():
            if v[2]==None:
                self._freezes_d[k]=(v[0], 'vary', False)
            else:
                self._links_d[k]=(v[0], 'expr', v[2])
        profile.disable()
        ps = pstats.Stats(profile)
        #ps.sort_stats('cumtime').print_stats(20)


    def _row_extract(self, id):
        """
        labels = self._labels_extract()
        ids = np.array([int(float(self._tab.GetCellValue(
                  i, np.where(labels == 'id')[0][0]))) \
                  for i in range(self._tab.GetNumberRows())])
        return np.where(id==ids)[0][0]
        """

        try:
            return np.where(id==self._ids)[0][0]
        except:
            return None


    def _text_colours(self):
        labels = self._labels_extract()
        for (r,c) in self._cells_sel:
            self._tab.SetCellTextColour(r, c, 'black')
        for i in range(self._tab.GetNumberRows()):
            idi = self._id_extract(i)
            #print(self._gui._sess_sel.systs)
            for m in self._gui._sess_sel.systs._mods_t:
                if idi in m['id']:
                    mod = m['mod']
                    #mod._pars.pretty_print()
                    #print(id(m))
            for p,v in mod._pars.items():
                if p.split('_')[-1] in ['z', 'logN', 'b', 'btur', 'resol']:
                    try:
                        c = np.where(labels==p.split('_')[-1])[0][0]
                    except:
                        c = None
                    r = i if c == 11 else self._row_extract(int(p.split('_')[-2]))
                    #if p.split('_')[-2] in ['45','46'] and p.split('_')[-1] == 'z':
                    #    print(idi, p,v)
                    if v.vary == False and r != None and c != None:
                        self._tab.SetCellTextColour(r, c, 'grey')
                    if v.expr != None and r != None and c != None:
                        r2 = self._row_extract(int(v.expr.split('_')[-2]))
                        c2 = np.where(labels==p.split('_')[-1])[0][0]
                        vs = v.expr.split('*')[0]
                        if vs not in self._links_c:
                            self._links_c[vs] = self._colours[self._colourc\
                                                    %len(self._colours)]
                            self._colourc += 1
                        self._tab.SetCellTextColour(r, c, self._links_c[vs])
                        self._tab.SetCellTextColour(r2, c2, self._links_c[vs])
            if not mod._active:
                cols = self._tab.GetNumberCols()
                for c in range(cols):
                    self._tab.SetCellTextColour(r, c, 'light grey')
