from .functions import get_selected_cells, trans_parse
from .gui_dialog import GUIDialogMiniSystems
from .vars import *
from collections import OrderedDict
import logging
import matplotlib.pyplot as plt
import numpy as np
import pprint
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
                    if n in ['logN', 'dlogN', 'b', 'db', 'resol', 'chi2r', \
                             'snr']:
                        format = '%3.3'
                    else:
                        format = '%3.7'
                    if np.abs(r[n])<1e-7 and r[n]!=0:
                        #self._tab.SetCellValue(j, i, "%3.7e" % r[n])
                        format += 'e'
                    else:
                        #self._tab.SetCellValue(j, i, "%3.7f" % r[n])
                        format += 'f'
                    self._tab.SetCellValue(j, i, format % r[n])
        self._tab.AutoSizeColumns(True)


    def _init(self, from_scratch=True):
        if not from_scratch:
            try:
                if self._tab.GetNumberRows() != 0:
                    self._tab.DeleteRows(0, self._tab.GetNumberRows())
            except:
                from_scratch = True
        if from_scratch:
            super(GUITable, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))

            self._panel = wx.Panel(self)
            self._tab = gridlib.Grid(self._panel)
            self._tab.CreateGrid(0, 0)
            self.SetPosition((0, wx.DisplaySize()[1]*0.5))

        coln = len(self._data.t.colnames)
        rown = len(self._data.t)-self._tab.GetNumberRows()
        self._tab.AppendCols(coln)
        self._tab.AppendRows(rown)


    def _init_data(self, from_scratch=True, autosort=True, attr=None):
        sess = self._gui._sess_sel
        if attr is None: attr = self._attr
        sess.json += self._gui._json_update("_tab", "_init_data",
                                            {"from_scratch": from_scratch,
                                             "autosort": autosort,
                                             "attr": attr})
        self._data = getattr(sess, attr)

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
        #self._gui._graph_det._graph._cursor_lines = []
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
        self._remove_data(row)

        sess.json += self._gui._json_update("cb", "_spec_update", {})
        sess.cb._spec_update()
        #self._tab.DeleteRows(pos=len(self._data.t), numRows=1)
        self._fill()
        self._gui._refresh(init_cursor=True)


    def _on_sort(self, event):
        labels = self._labels_extract()
        self._data.t.sort([labels[self._gui._col_sel], 'id'])
        self._gui._refresh(autosort=False)

    def _on_sort_reverse(self, event):
        labels = self._labels_extract()
        self._data.t['id'] = -1*self._data.t['id']
        self._data.t.sort([labels[self._gui._col_sel], 'id'], reverse=True)
        self._data.t['id'] = -1*self._data.t['id']
        self._gui._refresh(autosort=False)

    def _on_view(self, event=None, from_scratch=True, autosort=True):
        sess = self._gui._sess_sel
        self._init_data(from_scratch, autosort)
        if autosort:
            if 'z' in self._data.t.colnames: self._data.t.sort(['z','id'])
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
        sess = self._gui._sess_sel

        if self._attr == 'systs':
            sess.json += self._gui._json_update("cb", "_systs_remove",
                                                {"rem": [row]})
            sess.json += self._gui._json_update("cb", "_mods_recreate", {})

            sess.cb._systs_remove([row])
            sess.cb._mods_recreate()
        else:
            sess.json += self._gui._json_update("_tab", "_remove_data",
                                                {"row": row})
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


    def _on_close(self, event, **kwargs):
        super(GUITableLineList, self)._on_close(event, **kwargs)
        del self._gui._tab_lines._data
        self._menu.FindItemById(self._tab_id[1]).Check(False)
        #self._open['lines'] = False


    def _on_view(self, event, **kwargs):
        super(GUITableLineList, self)._on_view(event, **kwargs)
        #self._open['lines'] = True

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


    def _on_close(self, event, **kwargs):
        super(GUITableSpectrum, self)._on_close(event, **kwargs)
        del self._gui._tab_spec._data
        self._menu.FindItemById(self._tab_id[0]).Check(False)
        #self._open['spec'] = False


    def _on_view(self, event, **kwargs):
        super(GUITableSpectrum, self)._on_view(event, **kwargs)
        #self._open['spec'] = True


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


    def _id_extract(self, row):
        labels = self._labels_extract()
        return int(self._tab.GetCellValue(row, np.where(labels == 'id')[0][0]))


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
                if self._tab.GetCellTextColour(row, col) in self._colours:
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

    def _on_close(self, event, **kwargs):
        super(GUITableSystList, self)._on_close(event, **kwargs)
        del self._gui._tab_systs._data
        self._menu.FindItemById(self._tab_id[2]).Check(False)
        #self._open['systs'] = False


    def _on_detail(self, event, span=30):
        if event.GetRow() == -1: return
        if not hasattr(self._gui, '_graph_det'):
            from .gui_graph import GUIGraphDetail
            GUIGraphDetail(self._gui, init_ax=False)
        else:
            self._gui._graph_det._graph._fig.clear()

        row = self._data.t[event.GetRow()]
        z = row['z']
        series = trans_parse(row['series'])
        self._gui._graph_main._z_sel = z
        self._gui._graph_main._series_sel = series
        self._gui._graph_main._refresh(self._gui._sess_sel)
        self._gui._graph_det._update(series, z, hwin_def)
        if not hasattr(self._gui, '_dlg_mini_systems') \
            or self._gui._dlg_mini_systems == None:
            GUIDialogMiniSystems(self._gui, "System controls", series=row['series'], z=row['z'])
        else:
            self._gui._dlg_mini_systems._refresh(row['series'], row['z'])

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


    def _on_edit(self, event):
        row, col = event.GetRow(), event.GetCol()
        value = float(self._tab.GetCellValue(row, col))
        labels = self._labels_extract()
        self._data.t[labels[col]][row] = value
        id = self._id_extract(row)
        mod = self._mod_extract(row)
        try:
            vary = mod._pars['lines_voigt_%i_%s' % (id, labels[col])].vary
            expr = mod._pars['lines_voigt_%i_%s' % (id, labels[col])].expr
            mod._pars['lines_voigt_%i_%s' % (id, labels[col])].set(
                value=value, vary=vary, expr=expr)
        except:
            vary = mod._pars['psf_gauss_%i_%s' % (id, labels[col])].vary
            expr = mod._pars['psf_gauss_%i_%s' % (id, labels[col])].expr
            mod._pars['psf_gauss_%i_%s' % (id, labels[col])].set(
                value=value, vary=vary, expr=expr)


    def _on_fit(self, event):
        row = self._gui._tab_popup._event.GetRow()
        series = self._tab.GetCellValue(row, 1)
        z = float(self._tab.GetCellValue(row, 3))
        logN = float(self._tab.GetCellValue(row, 5))
        b = float(self._tab.GetCellValue(row, 7))
        cb = self._gui._sess_sel.cb
        mod = self._mod_extract(row)
        cb._syst_fit(mod, max_nfev_def)
        cb._spec_update()
        self._gui._refresh(init_cursor=True)


    def _on_freeze_par(self, event):
        popup = self._gui._tab_popup
        row = popup._event.GetRow()
        col = popup._event.GetCol()
        """
        for i in (33,34,35,40,41,42):
            try:
                print('Freeze before', self._freezes_d['lines_voigt_%i_b' % i])
            except:
                pass
        """
        for (r, c) in self._cells_sel:
            id, parn = self._key_extract(r, c)
            if self._tab.GetCellTextColour(row, col) != 'black':
                self._freezes_d[parn] = (id, 'vary', True)
            else:
                self._freezes_d[parn] = (id, 'vary', False)
        self._tab.ForceRefresh()
        systs = self._gui._sess_sel.systs
        """
        print('before constrain')
        for m in systs._mods_t:
            if (33 in m['id']):
                m['mod']._pars.pretty_print()
        """

        for v in self._freezes_d:
            if v in self._links_d and self._links_d[v][2] != '' and self._freezes_d[v][2] == True:
                #print(v, self._links_d[v], self._freezes_d[v])
                self._freezes_d[v] = (self._freezes_d[v][0],
                                      self._freezes_d[v][1], False)
                #print(v, self._links_d[v], self._freezes_d[v])

        systs._constrain(self._freezes_d)
        #systs._constrain(self._links_d)
        #self._gui._sess_sel.cb._mods_recreate2(only_constr=True)
        """
        print('after constrain')
        for m in systs._mods_t:
            if (33 in m['id']):
                m['mod']._pars.pretty_print()
        for i in (33,34,35,40,41,42):
            try:
                print('Freeze after', self._freezes_d['lines_voigt_%i_b' % i])
            except:
                pass
        """
        self._text_colours()


    def _on_freeze_par_all(self, event=None, col=None):
        if event is not None:
            col = self._gui._tab_popup._event.GetCol()
        for i in range(self._tab.GetNumberRows()):
            id, parn = self._key_extract(i, col)
            self._freezes_d[parn] = (id, 'vary', False)
        self._tab.ForceRefresh()
        self._gui._sess_sel.systs._constrain(self._freezes_d)
        self._text_colours()

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
            title = ['Sort ascending', 'Sort descending', 'sep', 'Histogram']
            attr = ['sort', 'sort_reverse', None, 'histogram']
            if col in [3,5,7]:
                title += ['sep', 'Freeze all']
                attr += [None, 'freeze_par_all']
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
        self._cells_sel = sorted(self._cells_sel, key=lambda tup: tup[0])
        ref = np.argmin([int(self._key_extract(r,c)[0]) for r,c in self._cells_sel])
        others = [self._cells_sel[c] for c in np.setdiff1d(range(len(self._cells_sel)), [ref])]
        id_r, val = self._key_extract(self._cells_sel[ref][0], self._cells_sel[ref][1])
        """
        for i in (33,34,35,40,41,42):
            try:
                print('Link before  ', self._links_d['lines_voigt_%i_b' % i])
            except:
                pass
        """
        for (r, c) in others:
            id, parn = self._key_extract(r, c)
            if self._tab.GetCellTextColour(r, c) in self._colours:#== 'forest green':
                self._links_d[parn] = (id, 'expr', '')
                self._links_d[val] = (id_r, 'expr', '')
            else:
                if parn != val:
                    self._tab.SetCellValue(r, c, self._tab.GetCellValue(self._cells_sel[ref][0], self._cells_sel[ref][1]))
                    self._links_d[parn] = (id, 'expr', val)
        self._tab.ForceRefresh()
        systs = self._gui._sess_sel.systs
        systs._constrain(self._links_d)
        self._gui._sess_sel.cb._mods_recreate2(only_constr=True)
        self._text_colours()

        """
        for m in systs._mods_t:
            if (33 in m['id']):
                m['mod']._pars.pretty_print()
        for i in (33,34,35,40,41,42):
            try:
                print('Link after  ', self._links_d['lines_voigt_%i_b' % i])
            except:
                pass
        """

    def _on_view(self, event, **kwargs):
        super(GUITableSystList, self)._on_view(event, **kwargs)
        #self._open['systs'] = True
        self._text_colours()
        self._tab.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK,
                       self._on_cell_right_click)
        for k, v in self._data._constr.items():
            if v[2]==None:
                self._freezes_d[k]=(v[0], 'vary', False)
            else:
                self._links_d[k]=(v[0], 'expr', v[2])


    def _row_extract(self, id):
        labels = self._labels_extract()
        ids = np.array([int(self._tab.GetCellValue(i, np.where(labels == 'id')[0][0])) \
               for i in range(self._tab.GetNumberRows())])
        return np.where(id==ids)[0][0]


    def _text_colours(self):
        labels = self._labels_extract()
        for (r,c) in self._cells_sel:
            self._tab.SetCellTextColour(r, c, 'black')
        for i in range(self._tab.GetNumberRows()):
            id = self._id_extract(i)
            for m in self._gui._sess_sel.systs._mods_t:
                if id in m['id']:
                    mod = m['mod']
            for p,v in mod._pars.items():
                if p.split('_')[-1] in ['z', 'logN', 'b', 'resol']:
                    c = np.where(labels==p.split('_')[-1])[0][0]
                    r = i if c == 9 else self._row_extract(int(p.split('_')[-2]))
                    if v.vary == False:
                        self._tab.SetCellTextColour(r, c, 'grey')
                    if v.expr != None:
                        r2 = self._row_extract(int(v.expr.split('_')[-2]))
                        c2 = np.where(labels==p.split('_')[-1])[0][0]
                        if v.expr not in self._links_c:
                            self._links_c[v.expr] = self._colours[self._colourc\
                                                    %len(self._colours)]
                            self._colourc += 1
                        self._tab.SetCellTextColour(r, c, self._links_c[v.expr])
                        self._tab.SetCellTextColour(r2, c2, self._links_c[v.expr])
