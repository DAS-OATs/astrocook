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
import time
import wx
import wx.grid as gridlib
import wx.lib.mixins.listctrl as listmix
import wx.lib.colourdb as cdb

"""
import scienceplots
from cycler import cycler
#plt.style.use('science')
#plt.rcParams.update({"font.family": "sans-serif", "font.size":14})
cc = cycler(plt.style.library['light']['axes.prop_cycle'])
cmap = [d['color'] for d in cc]
"""

max_rows = 2000

ds_x = int(wx.DisplaySize()[0]*0.98)
ds_y = int(wx.DisplaySize()[1]*0.88)
do_x = int(wx.DisplaySize()[0]*0.01)
do_y = int(wx.DisplaySize()[1]*0.05)

size_x_table = ds_x*3//5
size_y_table = ds_y//4
offset_x_table = do_x
offset_y_table = do_y+size_y_table


class GUITable(wx.Frame):
    """ Class for the GUI table frame """

    def __init__(self,
                 gui,
                 attr,
                 title="Table",
                 size_x=size_x_table,
                 size_y=size_y_table):

        self._gui = gui
        self._attr = attr
        self._title = title
        self._size_x = size_x
        self._size_y = size_y
        self._tab_id = self._gui._menu_tab_id
        self._gui._tab = self
        self._menu = self._gui._panel_sess._menu._view._menu
        super(GUITable, self).__init__(parent=None, title=self._title,
                                       size=(self._size_x, self._size_y))
        self._shown = False


    def _col_remove(self, cols, attr=None, log=True):
        sess = self._gui._sess_sel
        if attr is None: attr = self._attr

        tab = getattr(self._gui, '_tab_'+attr)
        coln = [self._data._t.colnames[c] for c in cols]
        tab._data.t.remove_columns(coln)


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
        if len(tab._data.t)>max_rows:
            logging.warning("The table is too long! I displayed only the first "
                            "{} rows. To display a different range, extract it "
                            "first.".format(max_rows))
            t = tab._data.t[:max_rows]
        else:
            t = tab._data.t
        for j, r in enumerate(t):
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
                        format += 'e'
                    else:
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
            tab.SetPosition((0, int(wx.DisplaySize()[1]*0.5)))

        coln = len(tab._data.t.colnames)
        rown = min(len(tab._data.t)-tab._tab.GetNumberRows(), max_rows)
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


    def _on_col_remove(self, event):
        cols = [self._gui._tab_popup._event.GetCol()]

        sess = self._gui._sess_sel
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
        sess.log.append_full('_tab', '_col_remove',
                             {'col': cols, 'attr': self._attr})
        sess.log.append_full('cb', '_spec_update', {})

        self._col_remove(cols, self._attr)
        self._gui._refresh(init_cursor=True)


    def _on_detail(self, event):
        if self._attr != 'systs': return
        if event.GetRow() == -1: return
        if not hasattr(self._gui, '_graph_det'):
            from .gui_graph import GUIGraphDetail
            GUIGraphDetail(self._gui, init_ax=False)
        elif len(self._gui._graph_det._graph._fig.axes) > 1:
            self._gui._graph_det._graph._fig.clear()
        size_x = int(wx.DisplaySize()[0]*0.4)
        size_y = int(wx.DisplaySize()[1]*0.4)
        self._gui._graph_det.SetSize(wx.Size(size_x, size_y))
        self._gui._graph_det._graph._init_ax(111)
        row = self._data.t[event.GetRow()]
        self._gui._sess_sel._xdet = row['x']
        self._gui._sess_sel._ydet = row['y']
        x = row['x']
        xlim, ylim = self._gui._graph_det._define_lim(x)
        self._gui._graph_split = False
        self._gui._graph_det._refresh(self._gui._sess_items, xlim=xlim,
                                      ylim=ylim)


    def _on_edit(self, event):
        row, col = event.GetRow(), event.GetCol()
        labels = self._labels_extract()
        label = labels[col]
        value = self._tab.GetCellValue(row, col)
        sess = self._gui._sess_sel
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
            GUIGraphHistogram(self._gui)
        else:
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
            title = ['Sort ascending', 'Sort descending', 'sep', 'Histogram', 'sep', 'Remove column']
            attr = ['sort', 'sort_reverse', None, 'histogram', None, 'col_remove']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title,
                                         attr), event.GetPosition())
        if col == -1:
            self.PopupMenu(GUITablePopup(self._gui, self, event, 'Remove row',
                                         'row_remove'), event.GetPosition())

    def _on_row_remove(self, event):
        if hasattr(self, '_cells_sel'):
            rows = [c[0] for c in self._cells_sel]
        else:
            rows = []
        if rows == []:
            rows = [self._gui._tab_popup._event.GetRow()]

        sess = self._gui._sess_sel
        if self._attr == 'systs':
            sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
            sess.log.append_full('cb', '_systs_remove', {'rem': rows})
        else:
            sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                        'from_scratch': False})
            sess.log.append_full('_tab', '_row_remove',
                                 {'row': rows, 'attr': self._attr})
        sess.log.append_full('cb', '_spec_update', {})

        self._row_remove(rows, self._attr)
        self._gui._refresh(init_cursor=True)


    def _on_replace(self, event):
        cb = self._gui._sess_sel.cb
        dlg = GUIDialogMethod(self._gui, 'Replace series', 'series_replace')
        self._gui._refresh(init_cursor=True)


    def _on_sort(self, event):
        labels = self._labels_extract()

        label = labels[self._gui._col_sel]
        sess = self._gui._sess_sel
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
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab', '_data_sort',
            {'label': label, 'attr': self._attr, 'reverse': True})
        self._data_sort(label=label, attr=self._attr, reverse=True)
        self._gui._refresh(autosort=False)


    def _on_view(self, event=None, from_scratch=True, autosort=False):
        self._view(event, from_scratch, autosort)


    def _row_remove(self, rows, attr=None, log=True):
        sess = self._gui._sess_sel
        if attr is None: attr = self._attr

        if self._attr == 'systs':
            rem_id = sess.cb._systs_remove(rows)
            sess.cb._spec_update()
        else:
            tab = getattr(self._gui, '_tab_'+attr)
            tab._data.t.remove_rows(rows)

    def _view(self, event=None, from_scratch=True, autosort=False):
        self._data_init(from_scratch, autosort)
        self._box = wx.BoxSizer(wx.VERTICAL)
        try:
            self._box.Add(self._tab, 1, wx.EXPAND)
        except:
            pass
        self._panel.SetSizer(self._box)
        self._tab.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self._on_edit)
        self._tab.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self._on_detail)
        self._tab.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK,
                       self._on_label_right_click)
        self.Bind(wx.EVT_CLOSE, self._on_close)
        self.Centre()
        self.SetPosition((offset_x_table, offset_y_table))
        self.Show()
        self._shown = True


class GUITableLineList(GUITable):
    """ Class for the GUI line list """

    def __init__(self,
                 gui,
                 title="Line table",
                 size_x=size_x_table,
                 size_y=size_y_table):

        super(GUITableLineList, self).__init__(gui, 'lines', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_lines = self


    def _on_close(self, event, **kwargs):
        super(GUITableLineList, self)._on_close(event, **kwargs)
        self._menu.FindItemById(self._tab_id[1]).Check(False)


    def _on_view(self, event, **kwargs):
        super(GUITableLineList, self)._on_view(event, **kwargs)


class GUITableModelList(GUITable):
    """ Class for the GUI model list """

    def __init__(self,
                 gui,
                 title="Model table",
                 size_x=size_x_table,
                 size_y=size_y_table):

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
                 size_x=size_x_table,
                 size_y=size_y_table):

        super(GUITableSpectrum, self).__init__(gui, 'spec', title, size_x,
                                               size_y)

        self._gui = gui
        self._gui._tab_spec = self


    def _on_close(self, event, **kwargs):
        super(GUITableSpectrum, self)._on_close(event, **kwargs)
        self._menu.FindItemById(self._tab_id[0]).Check(False)


    def _on_view(self, event, **kwargs):
        super(GUITableSpectrum, self)._on_view(event, **kwargs)


class GUITableSystList(GUITable):
    """ Class for the GUI system list """

    def __init__(self,
                 gui,
                 title="System table",
                 size_x=size_x_table,
                 size_y=size_y_table):

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

    def _col_idx_from_label(self, label_to_find):
        labels = self._labels_extract() # Assumes self._labels_extract() works
        try:
            return np.where(labels == label_to_find)[0][0]
        except IndexError:
            logging.error(f"Column with label '{label_to_find}' not found in table.")
            # Depending on usage, either raise an error or return a value like -1
            # For robust usage in _data_link_bt, it should probably raise or be checked.
            raise ValueError(f"Column '{label_to_find}' not found.")

    def _data_cell_right_click(self, row, col):
        sel = get_selected_cells(self._tab)
        if len(sel) == 1:
            self._tab.SetGridCursor(row, col)
            sel = get_selected_cells(self._tab)
        self._data_cells_desel()
        for s in sel:
            if s[1] in [1, 3, 5, 7, 9]:
                row = s[0]
                col = s[1]
                self._data_cells_sel(row, col)


    def _data_cells_sel(self, row, col, log=True):
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
        tt = time.time()
        if not hasattr(self._gui, '_graph_det'):
            from .gui_graph import GUIGraphDetail
            GUIGraphDetail(self._gui, init_ax=False)
        else:
            self._gui._graph_det._graph._fig.clear()
        z = row_z
        series = trans_parse(row_series)
        self._gui._graph_main._z_sel = z
        self._gui._graph_main._series_sel = series
        self._gui._graph_main._refresh(self._gui._sess_sel)
        self._gui._graph_det._update(series, z, hwin_def)
        if not hasattr(self._gui, '_dlg_mini_systems') \
            or self._gui._dlg_mini_systems == None:
            GUIDialogMiniSystems(self._gui, "System controls", series=row_series, z=row_z)
        else:
            try:
                self._gui._dlg_mini_systems._refresh(row_series, row_z)
            except:
                GUIDialogMiniSystems(self._gui, "System controls", series=row_series, z=row_z)
        dlg_mini_systems = self._gui._dlg_mini_systems
        dlg_mini_systems._menu.FindItemById(dlg_mini_systems._dlg_id[0]).Check(True)

        """
        # Color background of systems in the same group
        mods_sel = np.where([row_id in i \
                             for i in self._gui._sess_sel.systs._mods_t['id']])
        for j, r in enumerate(self._data.t):
            for i in range(len(self._data.t.colnames)):
                if r['id'] in np.array(self._gui._sess_sel.systs._mods_t['id'][mods_sel][0]):
                    self._tab.SetCellBackgroundColour(j, i, 'cyan')
                else:
                    self._tab.SetCellBackgroundColour(j, i, None)
        """
        s = self._gui._sess_sel.systs._d[row_id]
        for j, r in enumerate(self._data.t):
            if r['id']==row_id:
                self._tab.SetCellBackgroundColour(j, 3, 'spring green')
            elif r['id'] in np.array(s._group):
                self._tab.SetCellBackgroundColour(j, 3, 'yellow')
            else:
                self._tab.SetCellBackgroundColour(j, 3, None)

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


    def _data_freeze_par_old(self, row, col):
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
        systs._constrain(self._freezes_d, source="GUITableSystList._data_freeze_par")
        self._text_colours()


    def _data_freeze_par(self, row_clicked, col_clicked):
        """
        Toggles the freeze state of selected parameters.
        If a parameter is linked, freezing it removes the link.
        If a parameter is frozen, this unfreezes it.
        If a parameter is free, this freezes it.
        """
        changes_to_apply = {}
        systs = self._gui._sess_sel.systs

        logging.debug(f"GUITableSystList._data_freeze_par: ENTRY. Clicked_row={row_clicked}, col={col_clicked}. _cells_sel={self._cells_sel}")
        logging.debug(f"  Initial _freezes_d: {self._freezes_d}")
        logging.debug(f"  Initial _links_d: {self._links_d}")

        if not self._cells_sel:
            logging.warning("_data_freeze_par: No cells selected (_cells_sel is empty).")
            # Potentially use row_clicked, col_clicked if _cells_sel is empty after a single click
            # For now, assume _cells_sel is populated correctly by _data_cell_right_click
            if row_clicked != -1 and col_clicked != -1 : # Make sure it's a valid cell
                 try:
                    # Attempt to use the clicked cell if _cells_sel is empty
                    param_id_clicked, param_full_name_clicked = self._key_extract(row_clicked, col_clicked)
                    cells_to_process = [(row_clicked, col_clicked)]
                    logging.debug(f"  _cells_sel empty, using clicked cell: ({row_clicked},{col_clicked}) for {param_full_name_clicked}")
                 except Exception as e:
                    logging.error(f"  _cells_sel empty, and error extracting key from clicked cell ({row_clicked},{col_clicked}): {e}")
                    self._text_colours() # Refresh colors anyway
                    self._tab.ForceRefresh()
                    return
            else:
                self._text_colours()
                self._tab.ForceRefresh()
                return
        else:
            cells_to_process = self._cells_sel

        for (sel_r, sel_c) in cells_to_process:
            try:
                param_id, param_full_name = self._key_extract(sel_r, sel_c)
            except IndexError: # Happens if _key_extract fails (e.g. non-parameter cell)
                logging.warning(f"  Skipping cell ({sel_r},{sel_c}): Could not extract parameter key.")
                continue

            if self._tab.GetCellTextColour(sel_r, 0) == 'light grey': # System inactive
                logging.debug(f"  Skipping {param_full_name} at ({sel_r},{sel_c}): System inactive.")
                continue

            is_currently_frozen = param_full_name in self._freezes_d
            is_currently_linked = param_full_name in self._links_d

            if is_currently_frozen:
                # Action: UNFREEZE
                changes_to_apply[param_full_name] = (param_id, 'vary', True)
                del self._freezes_d[param_full_name]
                logging.info(f"  {param_full_name}: UNFREEZE. Removed from _freezes_d.")
            elif is_currently_linked:
                # Action: FREEZE (this overrides/removes the link)
                changes_to_apply[param_full_name] = (param_id, 'vary', False)
                self._freezes_d[param_full_name] = (param_id, 'vary', False) # Add to freezes
                del self._links_d[param_full_name] # Remove from links
                logging.info(f"  {param_full_name}: FREEZE (was linked). Added to _freezes_d, removed from _links_d.")
            else: # Currently free
                # Action: FREEZE
                changes_to_apply[param_full_name] = (param_id, 'vary', False)
                self._freezes_d[param_full_name] = (param_id, 'vary', False)
                logging.info(f"  {param_full_name}: FREEZE (was free). Added to _freezes_d.")

        if not changes_to_apply:
            logging.debug("  _data_freeze_par: No changes to apply.")
            self._text_colours()
            self._tab.ForceRefresh()
            return

        logging.debug(f"  _data_freeze_par: changes_to_apply: {changes_to_apply}")
        systs._constrain(changes_to_apply, source="GUITableSystList._data_freeze_par")

        logging.debug(f"  _data_freeze_par: Final _freezes_d: {self._freezes_d}")
        logging.debug(f"  _data_freeze_par: Final _links_d: {self._links_d}")

        self._text_colours()
        self._tab.ForceRefresh()
        logging.debug(f"GUITableSystList._data_freeze_par: EXIT")

    def _data_freeze_par_all_old(self, col, reverse):
        par = self._labels_extract()[col]
        self._tab.ForceRefresh()
        self._freezes_d = self._gui._sess_sel.systs._freeze_par(par, [], reverse)
        self._text_colours()
        if reverse:
            self._gui._refresh(init_cursor=True)


    def _data_freeze_par_all(self, col_to_process, freeze_action=True): # True to freeze, False to unfreeze
        """
        Freezes or unfreezes all applicable parameters in the specified column.
        """
        changes_to_apply = {}
        systs = self._gui._sess_sel.systs
        action_str = "FREEZE_ALL" if freeze_action else "UNFREEZE_ALL"
        
        logging.debug(f"GUITableSystList._data_freeze_par_all: ENTRY. Column_idx={col_to_process}, Action: {action_str}")
        logging.debug(f"  Initial _freezes_d: {self._freezes_d}")
        logging.debug(f"  Initial _links_d: {self._links_d}")
    
        if col_to_process < 0 or col_to_process >= self._tab.GetNumberCols():
            logging.error(f"  Invalid column index {col_to_process}.")
            return
    
        # param_name_in_col = self._labels_extract()[col_to_process] # e.g. 'z', 'logN', 'b'
    
        for r_idx in range(self._tab.GetNumberRows()):
            try:
                # We are iterating rows for a GIVEN column (col_to_process)
                param_id, param_full_name = self._key_extract(r_idx, col_to_process) 
            except IndexError:
                # This row/col might not be a valid parameter cell (e.g. header, or parse error)
                logging.warning(f"  Skipping cell ({r_idx},{col_to_process}) for {action_str}: Could not extract parameter key.")
                continue
                
            if self._tab.GetCellTextColour(r_idx, 0) == 'light grey': # System inactive
                logging.debug(f"  Skipping {param_full_name} for {action_str}: System inactive.")
                continue
            
            if freeze_action: # Intent is to FREEZE
                if param_full_name not in self._freezes_d: # Only if not already frozen
                    changes_to_apply[param_full_name] = (param_id, 'vary', False)
                    self._freezes_d[param_full_name] = (param_id, 'vary', False)
                    if param_full_name in self._links_d: # If it was linked, remove link
                        del self._links_d[param_full_name]
                    logging.info(f"  {param_full_name}: {action_str}. Added to _freezes_d, removed from _links_d if present.")
            else: # Intent is to UNFREEZE
                if param_full_name in self._freezes_d: # Only if currently frozen
                    changes_to_apply[param_full_name] = (param_id, 'vary', True)
                    del self._freezes_d[param_full_name]
                    logging.info(f"  {param_full_name}: {action_str}. Removed from _freezes_d.")
                # If unfreezing, we don't automatically re-link it if it was previously linked then frozen.
                # It just becomes free. User would have to re-link manually if desired.
    
        if not changes_to_apply:
            logging.debug(f"  _data_freeze_par_all: No changes to apply for {action_str}.")
            self._text_colours()
            self._tab.ForceRefresh()
            return
    
        logging.debug(f"  _data_freeze_par_all: changes_to_apply for {action_str}: {changes_to_apply}")
        systs._constrain(changes_to_apply, source=f"GUITableSystList._data_freeze_par_all ({action_str})")
    
        logging.debug(f"  _data_freeze_par_all: Final _freezes_d: {self._freezes_d}")
        logging.debug(f"  _data_freeze_par_all: Final _links_d: {self._links_d}")
        
        self._text_colours()
        self._tab.ForceRefresh()
        logging.debug(f"GUITableSystList._data_freeze_par_all: EXIT for {action_str}")


    def _data_init(self, from_scratch=True, autosort=False, attr=None):
        super(GUITableSystList, self)._data_init(from_scratch, autosort, attr)
        labels = self._labels_extract()
        self._ids = np.array([int(float(self._tab.GetCellValue(
                              i, np.where(labels == 'id')[0][0]))) \
                              for i in range(self._tab.GetNumberRows())])

    def _data_link_bt_old(self, row, col):
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
        systs._constrain(self._links_d, source="GUITableSystList._data_link_bt")
        self._gui._sess_sel.cb._mods_recreate2(only_constr=True)
        self._text_colours()


    def _data_link_bt(self, row_clicked, col_clicked):
        """
        Thermally links Doppler 'b' parameters of selected components (>=2) of the SAME ION
        to the 'b' parameter of the component with the smallest ID in the selection.
        Expression: b_other = b_ref * sqrt(mass_ion / mass_ion) = b_ref (if same ion)
        If different ions are selected, this should error or link b*sqrt(T) if that's the intent.
        The current name implies linking 'b' for thermal broadening, which requires temperature.
        If it's linking b values assuming same temperature: b_i ~ sqrt(1/m_i).
        b_other = b_ref * sqrt(m_ref / m_other)

        If already linked this way, unlinks them.
        If parameters are frozen, linking them unfreezes them.
        """
        changes_to_apply = {}
        systs = self._gui._sess_sel.systs

        logging.debug(f"GUITableSystList._data_link_bt: ENTRY. Clicked_row={row_clicked}, col={col_clicked}. _cells_sel={self._cells_sel}")
        logging.debug(f"  Initial _freezes_d: {self._freezes_d}")
        logging.debug(f"  Initial _links_d: {self._links_d}")

        if len(self._cells_sel) < 2:
            logging.warning("  _data_link_bt: Fewer than 2 cells selected.")
            wx.MessageBox("Select at least two 'b' parameter cells to link thermally.", "Linking Error", wx.OK | wx.ICON_WARNING)
            self._text_colours()
            self._tab.ForceRefresh()
            return

        selected_params_info = []
        param_col_label = self._labels_extract()[col_clicked] # Get label of clicked column, e.g. 'b'
        if param_col_label != 'b': # Ensure this is for 'b' parameters
            logging.error(f"  _data_link_bt: Clicked column '{param_col_label}' is not 'b'. This function is for 'b' params only.")
            wx.MessageBox(f"Thermal linking is for 'b' parameters only. Clicked on '{param_col_label}'.", "Linking Error", wx.OK | wx.ICON_ERROR)
            return

        for (sel_r, sel_c) in self._cells_sel:
            try:
                # Ensure we are only processing 'b' parameters from the selection
                if self._labels_extract()[sel_c] != 'b':
                    logging.warning(f"  Skipping cell ({sel_r},{sel_c}) for bT link: Not a 'b' parameter column.")
                    continue

                param_id, param_full_name = self._key_extract(sel_r, sel_c)
                if self._tab.GetCellTextColour(sel_r, 0) == 'light grey':
                    logging.debug(f"  Skipping {param_full_name} at ({sel_r},{sel_c}): System inactive.")
                    continue
                
                # Get series and ion for mass calculation
                series_str = self._tab.GetCellValue(sel_r, self._col_idx_from_label('series')) # Assuming 'series' column exists
                transitions = trans_parse(series_str)
                if not transitions:
                    logging.warning(f"  Skipping {param_full_name}: Could not parse transitions from series '{series_str}'.")
                    continue
                # The key for mass_d should be the resolved transition string itself,
                # assuming mass_d is built with these strings as keys.
                # We'll use the first resolved transition for this component's mass.
                mass_d_key = transitions[0] 

                if mass_d_key not in mass_d: # mass_d is from voigtastro.vars
                    logging.warning(f"  Skipping {param_full_name}: Resolved transition '{mass_d_key}' (from series '{series_str}') not found as a key in mass_d.")
                    continue
 
                ion_mass = mass_d[mass_d_key]
                # For logging or other purposes, you might still want a simplified ion name
                simple_ion_name = mass_d_key.split('_')[0] # Or a more robust extraction if needed
                
                selected_params_info.append({
                    'id': param_id, 'name': param_full_name, 'row': sel_r, 'col': sel_c,
                    'mass_d_key': mass_d_key, # Store the key used for mass_d
                    'ion_for_display': simple_ion_name, # For logging/display
                    'mass': ion_mass
                })
            except Exception as e: # Catch potential errors like _col_idx_from_label
                logging.error(f"  Error processing cell ({sel_r},{sel_c}) for bT link: {e}")
                continue

        if len(selected_params_info) < 2:
            logging.warning("  _data_link_bt: Fewer than 2 valid 'b' parameters from active systems selected.")
            wx.MessageBox("Select at least two valid 'b' parameter cells from active systems to link thermally.", "Linking Error", wx.OK | wx.ICON_WARNING)
            self._text_colours()
            self._tab.ForceRefresh()
            return

        selected_params_info.sort(key=lambda p: p['id'])
        ref_param = selected_params_info[0]
        other_params = selected_params_info[1:]

        logging.debug(f"  Reference param for bT linking: {ref_param['name']} (id: {ref_param['id']}, mass_key: {ref_param['mass_d_key']})")

        # Check for existing thermal link
        all_currently_linked_thermally = True
        if ref_param['name'] in self._links_d:
             all_currently_linked_thermally = False # Ref should not be dependent

        for p_other in other_params:
            mass_ratio_sqrt = np.sqrt(ref_param['mass'] / p_other['mass'])
            expected_expr = f"{ref_param['name']}*{mass_ratio_sqrt:.14f}" # Match precision from old code
            is_linked_as_expected = (
                p_other['name'] in self._links_d and
                self._links_d[p_other['name']][2] == expected_expr and
                self._links_d[p_other['name']][0] == p_other['id']
            )
            if not is_linked_as_expected:
                all_currently_linked_thermally = False
                break

        if all_currently_linked_thermally and len(other_params) > 0:
            logging.info(f"  Detected existing thermal link to {ref_param['name']}. Action: UNLINK all in selection.")
            for p_info in selected_params_info:
                if p_info['name'] in self._links_d:
                    changes_to_apply[p_info['name']] = (p_info['id'], 'expr', '')
                    del self._links_d[p_info['name']]
                    logging.info(f"    {p_info['name']}: UNLINK (was bT). Removed from _links_d.")
        else:
            logging.info(f"  Action: Thermally LINK other selected 'b' parameters to {ref_param['name']}.")
            if ref_param['name'] in self._links_d:
                changes_to_apply[ref_param['name']] = (ref_param['id'], 'expr', '')
                del self._links_d[ref_param['name']]
                logging.info(f"    Reference {ref_param['name']} was linked. Unlinking it.")
            if ref_param['name'] in self._freezes_d:
                changes_to_apply[ref_param['name']] = (ref_param['id'], 'vary', True)
                del self._freezes_d[ref_param['name']]
                logging.info(f"    Reference {ref_param['name']} was frozen. Unfreezing it.")

            for p_other in other_params:
                # The mass_ratio_sqrt is now calculated using the stored masses
                mass_ratio_sqrt = np.sqrt(ref_param['mass'] / p_other['mass'])
                link_expression = f"{ref_param['name']}*{mass_ratio_sqrt:.14f}"

                changes_to_apply[p_other['name']] = (p_other['id'], 'expr', link_expression)
                self._links_d[p_other['name']] = (p_other['id'], 'expr', link_expression)
                if p_other['name'] in self._freezes_d:
                    del self._freezes_d[p_other['name']]
                logging.info(f"    {p_other['name']} (mass_key {p_other['mass_d_key']}): LINK_BT to {link_expression}. Added to _links_d, removed from _freezes_d if present.")

                # No need for the p_other['ion'] == ref_param['ion'] check anymore if linking any two b's
                # based on their respective masses from mass_d, as the mass ratio handles it.
                # If you did want to restrict to same ion, you'd use 'ion_for_display' or re-extract.


        if not changes_to_apply:
            logging.debug("  _data_link_bt: No changes to apply.")
            self._text_colours()
            self._tab.ForceRefresh()
            return

        logging.debug(f"  _data_link_bt: changes_to_apply: {changes_to_apply}")
        systs._constrain(changes_to_apply, source="GUITableSystList._data_link_bt")

        logging.debug(f"  _data_link_bt: Final _freezes_d: {self._freezes_d}")
        logging.debug(f"  _data_link_bt: Final _links_d: {self._links_d}")

        self._text_colours()
        self._tab.ForceRefresh()
        logging.debug(f"GUITableSystList._data_link_bt: EXIT")

        # Helper to find column index (you might have this elsewhere)
        def _col_idx_from_label(self, label_to_find):
            labels = self._labels_extract()
            try:
                return np.where(labels == label_to_find)[0][0]
            except IndexError:
                logging.error(f"Column with label '{label_to_find}' not found.")
                raise # Or return -1 and handle

    def _data_link_par_old(self, row, col):
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
                    # Value from data table, to get full precision
                    vt = self._data.t[parn.split('_')[-1]]\
                             [self._cells_sel[ref][0]]
                    self._tab.SetCellValue(r, c, v)
                    self._links_d[parn] = (id, 'expr', val)
                    self._data_edit(r, labels[c], vt, update_mod=False)
        self._tab.ForceRefresh()
        systs = self._gui._sess_sel.systs
        systs._constrain(self._links_d, source="GUITableSystList._data_link_par")
        self._gui._sess_sel.cb._mods_recreate2(only_constr=True)
        self._text_colours()


    def _data_link_par(self, row_clicked, col_clicked):
        """
        Links selected parameters (>=2) to the one with the smallest ID in the selection.
        If parameters are frozen, linking them unfreezes them.
        If already linked (any param in selection to the ref), unlinks them.
        """
        changes_to_apply = {}
        systs = self._gui._sess_sel.systs

        logging.debug(f"GUITableSystList._data_link_par: ENTRY. Clicked_row={row_clicked}, col={col_clicked}. _cells_sel={self._cells_sel}")
        logging.debug(f"  Initial _freezes_d: {self._freezes_d}")
        logging.debug(f"  Initial _links_d: {self._links_d}")

        if len(self._cells_sel) < 2:
            logging.warning("  _data_link_par: Fewer than 2 cells selected. Linking requires at least two parameters.")
            wx.MessageBox("Select at least two parameter cells to link.", "Linking Error", wx.OK | wx.ICON_WARNING)
            self._text_colours()
            self._tab.ForceRefresh()
            return

        # Extract (param_id, param_full_name, row, col) for all selected cells
        selected_params_info = []
        for (sel_r, sel_c) in self._cells_sel:
            try:
                param_id, param_full_name = self._key_extract(sel_r, sel_c)
                if self._tab.GetCellTextColour(sel_r, 0) == 'light grey':
                    logging.debug(f"  Skipping {param_full_name} at ({sel_r},{sel_c}): System inactive.")
                    continue
                selected_params_info.append({'id': param_id, 'name': param_full_name, 'row': sel_r, 'col': sel_c})
            except IndexError:
                logging.warning(f"  Skipping cell ({sel_r},{sel_c}) for linking: Could not extract parameter key.")
                continue
            
        if len(selected_params_info) < 2:
            logging.warning("  _data_link_par: Fewer than 2 valid, active parameters selected.")
            wx.MessageBox("Select at least two valid parameter cells from active systems to link.", "Linking Error", wx.OK | wx.ICON_WARNING)
            self._text_colours()
            self._tab.ForceRefresh()
            return

        # Sort by param_id to find the reference parameter (smallest ID)
        selected_params_info.sort(key=lambda p: p['id'])
        ref_param = selected_params_info[0]
        other_params = selected_params_info[1:]

        logging.debug(f"  Reference param for linking: {ref_param['name']} (id: {ref_param['id']})")

        # Check if this set is ALREADY linked together (ref_param linked to itself, others linked to ref_param)
        # A simple check: is the reference parameter itself currently part of a link *as a dependent* or are others linked to it?
        # More robust: check if all 'other_params' are linked to 'ref_param' with the simple expression.
        all_currently_linked_to_ref = True
        if ref_param['name'] in self._links_d: # Ref param should not be a dependent if it's the reference
            all_currently_linked_to_ref = False

        for p_other in other_params:
            if not (p_other['name'] in self._links_d and \
                    self._links_d[p_other['name']][2] == ref_param['name'] and \
                    self._links_d[p_other['name']][0] == p_other['id']): # check id too
                all_currently_linked_to_ref = False
                break
            
        if all_currently_linked_to_ref and len(other_params) > 0 : # only consider unlinking if there were others
            # Action: UNLINK ALL in the selection
            logging.info(f"  Detected existing link to {ref_param['name']}. Action: UNLINK all in selection.")
            for p_info in selected_params_info: # Unlink ref as well if it was somehow linked
                if p_info['name'] in self._links_d:
                    changes_to_apply[p_info['name']] = (p_info['id'], 'expr', '') # Unlink
                    del self._links_d[p_info['name']]
                    logging.info(f"    {p_info['name']}: UNLINK. Removed from _links_d.")
        else:
            # Action: LINK OTHERS to REF_PARAM
            logging.info(f"  Action: LINK other selected parameters to {ref_param['name']}.")
            # Ensure ref_param is not linked to something else and not frozen
            if ref_param['name'] in self._links_d: # Ref param itself is a dependent! Cannot be.
                changes_to_apply[ref_param['name']] = (ref_param['id'], 'expr', '') # Unlink ref_param
                del self._links_d[ref_param['name']]
                logging.info(f"    Reference {ref_param['name']} was linked. Unlinking it.")
            if ref_param['name'] in self._freezes_d: # Ref param is frozen! Unfreeze it.
                changes_to_apply[ref_param['name']] = (ref_param['id'], 'vary', True) # Unfreeze
                del self._freezes_d[ref_param['name']]
                logging.info(f"    Reference {ref_param['name']} was frozen. Unfreezing it.")

            link_expression = ref_param['name'] # Simple link: P_other = P_ref

            for p_other in other_params:
                changes_to_apply[p_other['name']] = (p_other['id'], 'expr', link_expression)
                self._links_d[p_other['name']] = (p_other['id'], 'expr', link_expression) # Add to links
                if p_other['name'] in self._freezes_d: # If linked one was frozen, unfreeze it
                    # The 'expr' in changes_to_apply will tell _constrain to set vary=True
                    del self._freezes_d[p_other['name']] 
                logging.info(f"    {p_other['name']}: LINK to {link_expression}. Added to _links_d, removed from _freezes_d if present.")

                # Update cell value in table for immediate visual feedback (optional, _constrain + refresh should do it)
                # ref_val_str = self._tab.GetCellValue(ref_param['row'], ref_param['col'])
                # self._tab.SetCellValue(p_other['row'], p_other['col'], ref_val_str)


        if not changes_to_apply:
            logging.debug("  _data_link_par: No changes to apply.")
            self._text_colours()
            self._tab.ForceRefresh()
            return

        logging.debug(f"  _data_link_par: changes_to_apply: {changes_to_apply}")
        systs._constrain(changes_to_apply, source="GUITableSystList._data_link_par")

        logging.debug(f"  _data_link_par: Final _freezes_d: {self._freezes_d}")
        logging.debug(f"  _data_link_par: Final _links_d: {self._links_d}")

        self._text_colours()
        self._tab.ForceRefresh()
        logging.debug(f"GUITableSystList._data_link_par: EXIT")

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
        title = []
        attr = []
        if col in [3, 5, 7, 9]:
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
        if col in [1]:
            if len(self._cells_sel) > 1:
                title += ['Fit all systems...', 'Remove all']
                attr += ['syst_fit', 'row_remove']
            else:
                title += ['Fit system...', 'Fit group...', 'Extract group...',
                          'Remove row']
                attr += ['syst_fit', 'group_fit', 'group_extract', 'row_remove']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                           event.GetPosition())


    def _on_close(self, event, **kwargs):
        super(GUITableSystList, self)._on_close(event, **kwargs)
        self._menu.FindItemById(self._tab_id[2]).Check(False)


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
        if col<=1:
            value = str(self._tab.GetCellValue(row, col))
        else:
            value = float(self._tab.GetCellValue(row, col))
        sess = self._gui._sess_sel
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab_systs', '_data_edit',
                             {'row': row, 'label': label, 'value': value})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        self._data_edit(row, label, value)


    def _on_fit(self, event):
        # Deprecated
        row = self._gui._tab_popup._event.GetRow()
        sess = self._gui._sess_sel
        sess.log.append_full('_tab', '_data_init', {'attr': self._attr,
                                                    'from_scratch': False})
        sess.log.append_full('_tab_systs', '_data_fit', {'row': row})
        self._data_init(from_scratch=False, attr='systs')
        self._data_fit(row)
        self._gui._refresh(init_cursor=True)


    def _on_group_extract(self, event):
        row = self._gui._tab_popup._event.GetRow()
        id = self._id_extract(row)
        params = [{'id': id}]
        dlg = GUIDialogMethod(self._gui, 'Extract group', 'group_extract',
                              params_last=params)
        self._gui._refresh(init_cursor=True)


    def _on_group_fit(self, event):
        row = self._gui._tab_popup._event.GetRow()
        id = self._id_extract(row)
        params = [{'id': id, 'refit_n': 0, 'chi2rav_thres': 1e-2,
                   'max_nfev': max_nfev_def}]
        dlg = GUIDialogMethod(self._gui, 'Fit group', 'group_fit',
                              params_last=params)
        self._gui._refresh(init_cursor=True)


    def _on_syst_fit(self, event):
        rows = [c[0] for c in self._cells_sel]
        if rows == []:
            rows = [self._gui._tab_popup._event.GetRow()]
        ids = [self._id_extract(r) for r in rows]
        params = [{'ids': ids, 'refit_n': 0, 'chi2rav_thres': 1e-2,
                   'max_nfev': max_nfev_def}]
        dlg = GUIDialogMethod(self._gui, 'Fit system', 'syst_fit',
                              params_last = params)
        self._gui._refresh(init_cursor=True)


    def _on_systs_fit(self, event):
        params = [{'refit_n': 0, 'chi2rav_thres': 1e-2, 'max_nfev': max_nfev_def,
                   'sel_fit': False, '_mod': None}]
        dlg = GUIDialogMethod(self._gui, 'Fit systems', 'systs_fit',
                              params_last = params)
        self._gui._refresh(init_cursor=True)


    def _on_fit_dialog(self, event):
        # Deprecated
        row = self._data.t[self._gui._tab_popup._event.GetRow()]
        dlg = GUIDialogMethod(self._gui, 'Fit system', 'syst_fit',
                              params_last=[{'num': row._index+1}])

        self._gui._refresh(init_cursor=True)


    def _on_freeze_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        sess = self._gui._sess_sel
        sess.log.append_full('_tab_systs', '_data_freeze_par',
                             {'row': row, 'col': col})
        if hasattr(self._gui, '_dlg_mini_log') \
            and self._gui._dlg_mini_log._shown:
            self._gui._dlg_mini_log._refresh()
        logging.info("Calling _data_freeze_par FIRST TIME")
        self._data_freeze_par(row, col)
        logging.info("Returned from _data_freeze_par FIRST TIME")



    def _on_freeze_par_all(self, event=None, col=None, reverse=False):
        if event is not None:
            col = self._gui._tab_popup._event.GetCol()
            sess = self._gui._sess_sel
            sess.log.append_full('_tab_systs', '_data_freeze_par_all',
                                 {'col': col})
            if hasattr(self._gui, '_dlg_mini_log') \
                and self._gui._dlg_mini_log._shown:
                self._gui._dlg_mini_log._refresh()
            self._data_freeze_par_all(col, reverse)


    def _on_improve(self, event):
        dlg = GUIDialogMethod(self._gui, 'Improve systems', 'systs_improve')
        self._gui._refresh(init_cursor=True)


    def _on_label_right_click(self, event):
        row, col = event.GetRow(), event.GetCol()
        if row == -1 and col>0:
            self._data_top_label_right_click(col)
            title = ['Sort ascending', 'Sort descending']
            attr = ['sort', 'sort_reverse']
            if col==1:
                title += ['sep', 'Replace...']
                attr += [None, 'replace']
            if col > 1:
                title += ['sep', 'Histogram']
                attr += [None, 'histogram']
            if col in [3,5,7]:
                title += ['sep', 'Freeze all', 'Unfreeze all']
                attr += [None, 'freeze_par_all', 'unfreeze_par_all']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                           event.GetPosition())
        if row == -1 and col == -1:
            title = ['Fit all systems...']
            attr = ['systs_fit']
            self.PopupMenu(GUITablePopup(self._gui, self, event, title, attr),
                event.GetPosition())


    def _on_link_bt(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        sess = self._gui._sess_sel
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


    def _on_unfreeze_par_all(self, event=None, col=None):
        self._on_freeze_par_all(event, col, reverse=True)


    def _on_view(self, event, **kwargs):
        #profile = cProfile.Profile()
        #profile.enable()

        #super(GUITableSystList, self)._on_view(event, **kwargs)
        #self._text_colours()
        #self._tab.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK,
        #               self._on_cell_right_click)
        #for k, v in self._data._constr.items():
        #    if v[2]==None:
        #        self._freezes_d[k]=(v[0], 'vary', False)
        #    else:
        #        self._links_d[k]=(v[0], 'expr', v[2])

                # Calls _data_init, _fill, _init for the table grid itself
        super(GUITableSystList, self)._on_view(event, **kwargs) 

        # Initialize/clear _freezes_d and _links_d before repopulating.
        # This is CRITICAL to ensure they don't accumulate old states if _on_view is called multiple times.
        self._freezes_d = {}
        self._links_d = {}
        
        logging.debug(f"GUITableSystList._on_view: Initializing _freezes_d/_links_d from systs._constr.")
        if hasattr(self, '_data') and hasattr(self._data, '_constr'): # self._data is GUITable._attr ('systs'), so self._data is session.systs
            systs_constr = self._data._constr # This is session.systs._constr
            logging.debug(f"GUITableSystList._on_view: systs._constr to process: {systs_constr}")

            for k_constr, v_constr_tuple in systs_constr.items():
                # k_constr is the full parameter name (e.g., 'lines_voigt_X_z')
                # v_constr_tuple is (target_id, param_suffix, value_or_expr_string)
                
                if k_constr.endswith("_backup"): # Important: Skip any backup keys in _constr
                    logging.debug(f"GUITableSystList._on_view: Skipping backup key from _constr: {k_constr}")
                    continue

                # Ensure v_constr_tuple has the expected structure (3 elements)
                if not (isinstance(v_constr_tuple, (list, tuple)) and len(v_constr_tuple) == 3):
                    logging.warning(f"GUITableSystList._on_view: Malformed _constr entry for {k_constr}: {v_constr_tuple}. Skipping.")
                    continue
                
                target_id, param_suffix, val_expr_str = v_constr_tuple
                
                # Further validation: check if target_id and param_suffix seem valid for k_constr
                # (e.g., does k_constr contain target_id and param_suffix strings?)
                # This helps catch severely corrupted _constr entries.
                if not (str(target_id) in k_constr and param_suffix in k_constr):
                    logging.warning(f"GUITableSystList._on_view: Potential mismatch in _constr entry {k_constr} -> {v_constr_tuple}. Proceeding cautiously.")

                if val_expr_str is None: # Parameter is frozen in _constr
                    # _freezes_d expects: {'full_param_name': (id, 'vary', False)}
                    self._freezes_d[k_constr] = (target_id, 'vary', False)
                elif val_expr_str == '': 
                    # Parameter was explicitly unlinked and is just varying; not actively constrained by link/freeze.
                    # It should NOT be in _freezes_d or _links_d.
                    # If it's in _constr with ('id', 'suffix', ''), it means _constrain put it there.
                    # This state means "actively set to free after being constrained".
                    # We usually represent "free" by absence from _freezes_d and _links_d.
                    pass # Do not add to _freezes_d or _links_d
                else: # Parameter has an expression string, so it's linked in _constr
                    # _links_d expects: {'full_param_name': (id, 'expr', 'expression_string')}
                    self._links_d[k_constr] = (target_id, 'expr', val_expr_str)
            
            logging.debug(f"GUITableSystList._on_view: Initialized _freezes_d: {self._freezes_d}")
            logging.debug(f"GUITableSystList._on_view: Initialized _links_d: {self._links_d}")
        else:
            logging.warning("GUITableSystList._on_view: self._data or self._data._constr not found. _freezes_d/_links_d may be empty.")

        # _text_colours() should be called AFTER _freezes_d and _links_d are populated,
        # if it relies on them. However, _text_colours primarily reads from lmfit model._pars.
        # The lmfit model._pars should already be in sync with _systs_constr due to prior operations
        # or model loading/recreation.
        # The order here is okay: super()._on_view sets up the table cells, then we sync our
        # GUI constraint dicts, then _text_colours applies colors based on lmfit state.
        self._text_colours()

        # Event binding is fine here.
        self._tab.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self._on_cell_right_click)
        
        # The old loop `for k, v in self._data._constr.items():` was almost correct.
        # The main additions are:
        # 1. Clearing _freezes_d and _links_d at the start.
        # 2. Skipping "_backup" keys.
        # 3. Handling the `val_expr_str == ''` case (explicitly free/unlinked) by not adding to either dict.
        # 4. Using the correct tuple structure for _freezes_d and _links_d.

        #profile.disable()
        #ps = pstats.Stats(profile)


    def _row_extract(self, id):
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
            for m in self._gui._sess_sel.systs._mods_t:
                if idi in m['id']:
                    mod = m['mod']
            for p,v in mod._pars.items():
                if p.split('_')[-1] in ['z', 'logN', 'b', 'btur', 'resol']:
                    try:
                        c = np.where(labels==p.split('_')[-1])[0][0]
                    except:
                        c = None
                    r = i if c == 11 else self._row_extract(int(p.split('_')[-2]))
                    if r!=None and c!=None:

                        if v.vary==False: # and r != None and c != None:
                            self._tab.SetCellTextColour(r, c, 'grey')
                        if v.expr!=None: #and r != None and c != None:
                            r2 = self._row_extract(int(v.expr.split('_')[-2]))
                            c2 = np.where(labels==p.split('_')[-1])[0][0]
                            vs = v.expr.split('*')[0]
                            if vs not in self._links_c:
                                self._links_c[vs] = self._colours[self._colourc\
                                                        %len(self._colours)]
                                self._colourc += 1
                            self._tab.SetCellTextColour(r, c, self._links_c[vs])
                            if r2 != None and c2 != None:
                                self._tab.SetCellTextColour(r2, c2, self._links_c[vs])
                        #if v.vary==True and v.expr==None:
                        #    self._tab.SetCellTextColour(r, c, 'black')
            if not mod._active:
                cols = self._tab.GetNumberCols()
                for c in range(cols):
                    self._tab.SetCellTextColour(r, c, 'light grey')
