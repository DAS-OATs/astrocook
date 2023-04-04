from . import *
from .defaults import *
from .functions import expr_eval, x_convert, get_home_limits, update_home_limits
from .gui_graph import *
from .gui_image import *
from .gui_log import *
from .gui_menu import *
from .gui_table import *
from .message import *
from astropy import table as at
from astropy import constants as ac
from collections import OrderedDict
from copy import deepcopy as dc
import datetime as dt
import json
import logging
from matplotlib import pyplot as plt
import numpy as np
from sphinx.util import docstrings as ds
import wx
import wx.lib.mixins.listctrl as listmix

class GUI(object):
    """ Class for the GUI. """

    def __init__(self, paths=None, flags=None):
        """ Constructor """

        self._flags = flags
        try:
            banner = 'ASTROCOOK 🍪 v'
            l = ['─']*(3+len(banner)+len(version))
            print("┌%s┐" % ''.join(l))
            print("│ %s%3s │" % (banner, version))
            print("└%s┘" % ''.join(l))
        except:
            l = ['-']*(17+len(version))
            print(''.join(l))
            print(" ASTROCOOK  v%3s " % version)
            print(''.join(l))
        print("Cupani et al. 2017-%s * INAF-OATs" % current_year)
        self._sess_list = []
        self._sess_item_list = []
        self._sess_sel = None
        self._sess_item_sel = []
        self._menu_spec_id = []
        self._menu_y_conv_id = []
        self._menu_lines_id = []
        self._menu_cont_id = []
        self._menu_nodes_id = []
        self._menu_systs_id = []
        self._menu_z0_id = []
        self._menu_mods_id = []
        self._menu_feats_id = []
        self._menu_tab_id = []
        self._defs = Defaults(self)
        self._panel_sess = GUIPanelSession(self)
        self._id_zoom = 9
        self._data_lim = None
        self._tag = ""
        self._ok = True
        GUIGraphMain(self)
        GUITableSpectrum(self)
        GUITableLineList(self)
        GUITableSystList(self)
        GUITableModelList(self)
        GUIImageCompleteness(self)
        GUIImageCorrectness(self)
        if paths == None:
            logging.info("Welcome! Try Session > Open...")
        else:
            logging.info("Welcome!")
            for p in paths:
                self._panel_sess._open_path = p
                if p[-4:] == 'json':
                    self._panel_sess._open_rec = 'json_load'
                    self._panel_sess.json_load(os.path.realpath(p))
                else:
                    self._panel_sess._open_rec = '_on_open'
                    self._panel_sess._on_open(os.path.realpath(p))

    def _flags_cond(self, flag):
        return self._flags is not None \
            and flag in [f[:len(flag)] for f in self._flags]

    def _flags_extr(self, flag):
        extr = [f.split('=')[-1] for f in self._flags if f[:len(flag)]==flag]
        if len(extr)>1:
            logging.warning("You gave me too many %s flags! I will consider "\
                            "only the first one." % flag)
        return extr[0]


    def _log_rerun(self, log):
        i = self._sess_list.index(self._sess_sel)
        self._panel_sess._tab.DeleteItem(i)
        del self._sess_list[i]
        del self._sess_item_list[i]

        # Run selected log
        self._log_run(json.loads(log), skip_tab=True)


    def _log_run(self, load, skip_tab=False):

        if hasattr(self, '_graph_det'):
            self._graph_det._on_close()

        for r in load['set_menu']:

            if r['cookbook']=='_dlg_mini_graph' and 'value' in r['params']:
                i = self._sess_list.index(self._sess_sel)
                r['params']['value'] = \
                    ('%i,'%i).join(r['params']['value'].split('SESS_SEL,'))

            if r['cookbook'][:8]=='cookbook' or r['cookbook']=='cb':
                cb = self._sess_sel.cb
            elif r['cookbook'] == '':
                cb = self
            else:
                rs = r['cookbook'].split('.')
                cb = getattr(self, rs[0])
                for s in rs[1:]:
                    cb = getattr(cb, s)

            skip = False
            if hasattr(cb, '_menu_view'):
                try:
                    if getattr(self, '_tab_'+r['params']['obj']+'_shown') \
                        and r['params']['check']==True:
                        skip = True
                except:
                    pass

            if not skip:
                if r['recipe'][0] != '_':
                    logging.info("I'm launching %s..." % r['recipe'])
                start = dt.datetime.now()
                out = getattr(cb, r['recipe'])(**r['params'])
                end = dt.datetime.now()
                if r['recipe'][0] != '_':
                    logging.info("I completed %s in %3.3f seconds!" \
                                 % (r['recipe'], (end-start).total_seconds()))
            if out is not None and out != 0:
                self._panel_sess._on_add(out, open=False)

            if out is None or out==0:
                if hasattr(cb, '_tag') and cb._tag=='cb':
                    self._sess_sel.log.append_full(cb._tag, r['recipe'],
                                                   r['params'])
                if hasattr(cb, '_tag') and cb._tag=='_panel_sess' \
                    and r['recipe']=='equalize':
                    sess_list = [self._sess_list[s] \
                                 for s in r['params']['_sel']]
                    self._sess_sel.log.merge_full(cb._tag, r['recipe'],
                                                  r['params'], sess_list,
                                                  self._sess_sel)
                if hasattr(cb, '_tab_id'):
                    try:
                        attr = r['params']['attr']
                    except:
                        attr = cb._attr
                    tab = '_tab_'+attr
                    if cb._attr == 'systs':
                        tab_tag = '_tab_systs'
                    else:
                        tab_tag = '_tab'
                    self._sess_sel.log.append_full(tab_tag, r['recipe'],
                                                   r['params'])
                if cb!=self and hasattr(cb, '_menu_view'):
                    self._sess_sel.log.append_full('_menu_view', r['recipe'],
                                                   r['params'])
            else:
                if hasattr(cb, '_tag') and cb._tag=='cb':
                    sess_list = [self._sess_list[sel_old]]
                    self._sess_sel.log.merge_full(cb._tag, r['recipe'],
                                                  r['params'], sess_list,
                                                  self._sess_sel)
                if hasattr(cb, '_tag') and cb._tag=='_panel_sess' \
                    and r['recipe']=='combine':
                    sess_list = [self._sess_list[s] \
                                 for s in r['params']['_sel']]
                    self._sess_sel.log.merge_full(cb._tag, r['recipe'],
                                                  r['params'], sess_list,
                                                  self._sess_sel)


            sel_old = self._sess_list.index(self._sess_sel)


    def _refresh(self, init_cursor=False, init_tab=True, init_bar=False,
                 autolim=True, autosort=True, _xlim=None, _ylim=None):
        """ Refresh the GUI after an action """
        self._defs = self._sess_sel.defs
        try:
            old_xaxis_unit = self._sess_sel._gui._graph_main._graph._xunit
        except:
            old_xaxis_unit = au.nm

        self._panel_sess._refresh()
        self._panel_sess._menu._refresh(init_bar=init_bar)
        if hasattr(self, '_dlg_mini_defs') \
            and self._dlg_mini_defs._shown:
            self._dlg_mini_defs._refresh()
        if hasattr(self, '_dlg_mini_graph') \
            and self._dlg_mini_graph._shown:
            self._graph_main._elem = self._sess_sel._graph_elem
            if hasattr(self, '_graph_det') \
                and hasattr(self._graph_det, '_elem'):
                self._graph_det._elem = self._sess_sel._graph_elem
            self._graph_main._lim = self._sess_sel._graph_lim
            self._dlg_mini_graph._refresh()
        else:
            if hasattr(self._sess_sel, '_graph_elem'):
                self._graph_main._elem = self._sess_sel._graph_elem
                if hasattr(self, '_graph_det') \
                    and hasattr(self._graph_det, '_elem'):
                    self._graph_det._elem = self._sess_sel._graph_elem
            else:
                self._graph_main._elem = elem_expand(graph_elem,
                    self._panel_sess._sel)
                if hasattr(self, '_graph_det') \
                    and hasattr(self._graph_det, '_elem'):
                    self._graph_det._elem = elem_expand(graph_elem,
                        self._panel_sess._sel)
            if hasattr(self._sess_sel, '_graph_lim'):
                self._graph_main._lim = self._sess_sel._graph_lim
            else:
                self._graph_main._lim = graph_lim_def

        if hasattr(self, '_dlg_mini_log') and self._dlg_mini_log._shown:
            self._dlg_mini_log._refresh()
        if hasattr(self, '_dlg_mini_meta') and self._dlg_mini_meta._shown:
            self._dlg_mini_meta._refresh()

        ax = self._graph_main._graph._ax
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        dl = self._data_lim
        if dl is None or ((xlim[0] <= dl[0] or dl[0]==0.) \
                      and (xlim[1] >= dl[1] or dl[1]==1.) \
                      and (ylim[0] <= dl[2] or dl[2]==0.) \
                      and (ylim[1] >= dl[3] or dl[3]==1.)):
            self._graph_main._graph._zoom = False
            dl = (xlim[0], xlim[1], ylim[0], ylim[1])

        # Convert xlim and dl if needed
        if not hasattr(self._sess_sel.spec, '_xunit_old'):
            self._sess_sel.spec._xunit_old = au.nm
            self._sess_sel.spec._zem = 0
        try:
            xunit_old = self._sess_sel.spec._xunit_old
            xunit = self._sess_sel.spec._xunit
            zem = self._sess_sel.spec._zem
            if xunit_old != xunit:
                xlim = (x_convert((xlim[0]*xunit_old), zem, xunit).value,
                        x_convert(xlim[1]*xunit_old, zem, xunit).value)
                dl = (x_convert(dl[0]*xunit_old, zem, xunit).value,
                      x_convert(dl[1]*xunit_old, zem, xunit).value, dl[2], dl[3])
        except:
            pass




        goodlim = True
        if xlim == (0.0, 1.0) and ylim == (0.0, 1.0):
            goodlim = False
        if autolim and goodlim and _xlim != None and _ylim != None:
            self._graph_main._refresh(self._sess_items, xlim=list(_xlim),
                                      ylim=list(_ylim))
        elif autolim and goodlim and _xlim != None:
            self._graph_main._refresh(self._sess_items, xlim=list(_xlim))
        elif autolim and goodlim and self._graph_main._graph._zoom:
            self._graph_main._refresh(self._sess_items, xlim=xlim, ylim=ylim)
        else:
            self._graph_main._refresh(self._sess_items)

        # update stack to fix bug with home button
        # get stack
        stack_elems = self._sess_sel._gui._graph_main._graph._toolbar._nav_stack._elements
        xaxis_unit = self._sess_sel._gui._graph_main._graph._xunit
        # get limits
        main_graph_home_limits = get_home_limits(stack_elems)
        if main_graph_home_limits is not None:
            mghl_x = main_graph_home_limits[0:2]
            # technically we should do the same for y
            mghl_y = main_graph_home_limits[2:]
            try:
                mghl_x = (x_convert(mghl_x[0]*old_xaxis_unit, zem, xaxis_unit).value,
                            x_convert(mghl_x[1]*old_xaxis_unit, zem, xaxis_unit).value)
            except:
                pass

    # update firts element of stack
    update_home_limits(stack_elems, xlim=mghl_x, ylim=None)

        if hasattr(self, '_graph_det'):
            if self._sess_sel.systs is None:
                self._graph_det._on_close()
            else:
                graph = self._graph_det._graph
                if hasattr(graph, '_axes'):
                    for key in graph._zems:
                        xunit = self._sess_sel.spec.x.unit
                        self._sess_sel.cb.x_convert(zem=graph._zems[key])
                        graph._ax = graph._axes[key]
                        xlim_det = graph._ax.get_xlim()
                        ylim_det = graph._ax.get_ylim()
                        if autolim or True:
                            self._graph_det._refresh(self._sess_items, text=key,
                                                     xlim=xlim_det, ylim=ylim_det,
                                                     init_cursor=init_cursor)
                        else:
                            self._graph_det._refresh(self._sess_items, text=key,
                                                     init_cursor=init_cursor)
                        init_cursor = False
                        self._sess_sel.cb.x_convert(zem=graph._zems[key], xunit=xunit)
                else:
                    xlim_det = graph._ax.get_xlim()
                    ylim_det = graph._ax.get_ylim()
                    if autolim:
                        self._graph_det._refresh(self._sess_items, xlim=xlim_det,
                                                 ylim=ylim_det,
                                                 init_cursor=init_cursor)
                    else:
                        self._graph_det._refresh(self._sess_items,
                                                 init_cursor=init_cursor)
        for s in ['spec', 'lines', 'systs']:
            if hasattr(self, '_tab_'+s) and init_tab:
                if hasattr(getattr(self, '_tab_'+s), '_data'):
                    if hasattr(getattr(self._sess_sel, s), '_t'):
                        index = ['spec', 'lines', 'systs'].index(s)
                        item = self._menu_view._menu.FindItemById(self._menu_tab_id[index])
                        view = item.IsChecked()
                        if view:
                            getattr(self, '_tab_'+s)._view(
                                event=None, from_scratch=False, autosort=autosort)
                    else:
                        try:
                            getattr(self, '_tab_'+s).Destroy()
                        except:
                            pass

                if hasattr(self, '_col_sel') \
                    and self._col_sel < self._col_tab.GetNumberCols():
                    try:
                        self._col_values = \
                            [float(self._col_tab.GetCellValue(i, self._col_sel)) \
                            for i in range(self._col_tab.GetNumberRows())]
                    except:
                        self._col_values = \
                            [self._col_tab.GetCellValue(i, self._col_sel) \
                            for i in range(self._col_tab.GetNumberRows())]
        if hasattr(self, '_tab_systs') and self._tab_systs._shown:
            try:
                self._tab_systs._text_colours()
            except:
                pass
        if hasattr(self, '_graph_hist'):
            self._graph_hist._refresh(self._sess_items)


    def _refresh_graph_det(self, init_cursor=False, autolim=True):
        graph = self._graph_det._graph
        if hasattr(graph, '_axes'):
            for key in graph._zems:
                xunit = self._sess_sel.spec.x.unit
                self._sess_sel.cb.x_convert(zem=graph._zems[key])
                graph._ax = graph._axes[key]
                xlim_det = graph._ax.get_xlim()
                ylim_det = graph._ax.get_ylim()
                if autolim:
                    self._graph_det._refresh(self._sess_items, text=key,
                                             xlim=xlim_det, ylim=ylim_det,
                                             init_cursor=init_cursor)
                else:
                    self._graph_det._refresh(self._sess_items, text=key,
                                             init_cursor=init_cursor)
                init_cursor = False
                self._sess_sel.cb.x_convert(zem=graph._zems[key], xunit=xunit)
        else:
            xlim_det = graph._ax.get_xlim()
            ylim_det = graph._ax.get_ylim()
            if autolim:
                self._graph_det._refresh(self._sess_items, xlim=xlim_det,
                                         ylim=ylim_det,
                                         init_cursor=init_cursor)
            else:
                self._graph_det._refresh(self._sess_items,
                                         init_cursor=init_cursor)

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

    def _get_selected_items(self):
        sel = []

        # start at -1 to get the first selected item
        current = -1
        while True:
            next = self.GetNextItem(current, wx.LIST_NEXT_ALL,
                                    wx.LIST_STATE_SELECTED)
            if next == -1:
                return sel

            sel.append(next)
            current = next

    def _insert_string_item(self, *args):
        self.InsertItem(*args)
        listmix.TextEditMixin.__init__(self)


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
        self.SetPosition((wx.DisplaySize()[0]*0.02, wx.DisplaySize()[0]*0.02))

        self._tag = "_panel_sess"

        # Import GUI
        self._gui = gui
        self._gui._panel_sess = self

        # Create table
        panel = wx.Panel(self)
        self._tab = GUIControlList(panel, 0)
        self._tab.InsertColumn(0, "name", width=330)
        self._tab.InsertColumn(1, "id", width=30)
        self._tab.InsertColumn(2, "active range", width=200)
        self._tab.InsertColumn(3, "# rows", width=100)
        self._tab.InsertColumn(4, "# lines", width=100)
        self._tab.InsertColumn(5, "# systems", width=100)
        self._tab.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self._on_veto)
        self._tab.Bind(wx.EVT_LIST_END_LABEL_EDIT, self._on_edit)
        self._tab.Bind(wx.EVT_LIST_ITEM_SELECTED, self._on_select)
        self._tab.Bind(wx.EVT_LIST_ITEM_DESELECTED, self._on_deselect)
        self._tab.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self._on_right_click)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        panel.SetSizer(self._box)
        self._menu = GUIMenu(self._gui)
        self.SetMenuBar(self._menu.bar())
        self.Show()
        self.Bind(wx.EVT_CLOSE, self._on_close)


    def _on_add(self, sess, open=True):
        # _sel is the last selection; _items is the list of all selections.
        sess.log = GUILog(self._gui)
        sess.defs = Defaults(self._gui)
        self._gui._defs = sess.defs

        missing = []
        for i in range(self._tab.GetItemCount()+1):
            if i not in self._gui._sess_item_list:
                missing.append(i)

        self._sel = missing[0]
        self._items = [self._sel]
        self._tab._insert_string_item(self._sel, "%s" % sess.name)
        self._tab.SetItem(self._tab.GetItemCount()-1, 1, "%s" % str(self._sel))
        self._gui._sess_list.append(sess)
        self._gui._sess_item_list.append(self._sel)

        sel_sort = np.argsort(self._gui._sess_item_list)
        self._gui._sess_list = list(np.array(self._gui._sess_list)[sel_sort])
        self._gui._sess_item_list = \
            list(np.array(self._gui._sess_item_list)[sel_sort])

        # Similarly, _sess_sel contains the last selected session; _sess_items
        # contains all selected sessions
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_items = [self._gui._sess_sel]
        if open:
            ko = self._gui._sess_sel.open()
        else:
            ko = False

        if not ko:
            x = sess.spec._safe(sess.spec.x)
            self._gui._sess_sel._graph_elem = elem_expand(graph_elem, self._sel)
            self._gui._sess_sel._graph_lim = graph_lim_def

            self._gui._refresh(init_tab=False, autolim=False)

            # Enable import from depending on how many sessions are present
            for menu in [self._menu._edit, self._menu._cb_general]:
                menu_dict = menu._menu.__dict__
                for m in menu_dict:
                    menu._menu.Enable(menu_dict[m]['start_id'],
                    menu._enable(menu_dict[m]['func'], menu_dict[m]['value']))
        self._gui._ok = not ko


    def _on_edit(self, event):
        self._gui._sess_list[self._sel].spec.meta['object'] = event.GetLabel()


    def _on_open(self, path, _flags=None):
        """ Behaviour for Session > Open """

        if _flags is None or _flags==[] and self._gui._flags is not None:
            _flags = self._gui._flags
        elif _flags is not None:
            self._gui._flags = _flags

        name = '.'.join(path.split('/')[-1].split('.')[:-1])

        if not os.path.exists(path):
            logging.warning("File %s doesn't exist. I'm ignoring it." % path)
            return 0

        logging.info("I'm loading file %s into a new session..." % path)
        sess = Session(gui=self._gui, path=path, name=name)
        self._gui._panel_sess._on_add(sess, open=True)
        sess.log.append_full('_panel_sess', '_on_open',
                             {'path': path, '_flags': _flags})

        if sess._open_twin:
            logging.info("I'm loading twin session %s..." % path)
            sess = Session(gui=self._gui, path=path, name=name, twin=True)
            self._gui._panel_sess._on_add(sess, open=True)

        if _flags is not None and '-s' in _flags:
            logging.info("I'm loading session for slice 0...")
            sess = Session(gui=self._gui, path=path, name='%s_0' \
                           % name, slice=0)
            self._gui._panel_sess._on_add(sess, open=True)
            logging.info("I'm loading session for slice 1...")
            sess = Session(gui=self._gui, path=path, name='%s_1' \
                           % name, slice=1)
            self._gui._panel_sess._on_add(sess, open=True)

        if _flags is not None and '-o' in [f[:2] for f in _flags]:

            # Select orders based on flag
            flag = np.array(_flags)[np.where(['-o' in f[:2] for f in _flags])][0][3:]
            if flag != '':
                olim = [int(f) for f in flag.split(',')]
                rlim = [np.where(sess._order==olim[0])[0][0],
                    np.where(sess._order==olim[1])[0][-1]]
                sess._row = rlim[0]
                lim = rlim[1]
            else:
                lim = np.inf

            # Create new sessions for the selected orders
            while sess._row is not None and sess._row <= lim:
                logging.info("I'm loading session for order %i, slice 0..." \
                    % sess._order[sess._row])
                sess = Session(gui=self._gui, path=path, name='%s_%i-0' \
                               % (name, sess._order[sess._row]), row=sess._row)
                self._gui._panel_sess._on_add(sess, open=True)
                logging.info("I'm loading session for order %i, slice 1..." \
                             % sess._order[sess._row])
                sess = Session(gui=self._gui, path=path, name='%s_%i-1' \
                               % (name, sess._order[sess._row]), row=sess._row)
                self._gui._panel_sess._on_add(sess, open=True)


    def _entry_select(self):
        try:
            self._gui._sess_sel.cb.sess = self._gui._sess_sel
        except:
            pass

        # Enable session equalize/combine depending on how many sessions are selected
        for menu in [self._menu._edit, self._menu._cb_general]:
            menu_dict = menu._menu.__dict__
            for m in menu_dict:
                menu._menu.Enable(menu_dict[m]['start_id'],
                menu._enable(menu_dict[m]['func'], menu_dict[m]['value']))

        item = self._tab.GetFirstSelected()

        # Selection via JSON
        if item == -1:
            self._gui._refresh()

        # Manual selection
        else:
            self._items = []
            while item != -1:
                self._items.append(item)
                item = self._tab.GetNextSelected(item)
                self._gui._sess_items = [self._gui._sess_list[i] for i in self._items]
            if self._gui._sess_item_sel != []:
                self._gui._refresh()


    def _on_close(self, event):
        logging.info("Bye!")
        self.Destroy()
        os._exit(1)


    def _on_deselect(self, event):
        self._sel = event.GetIndex()
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_item_sel = []


    def _on_rerun(self, event):
        from .gui_dialog import GUIDialogMiniLog
        self._gui._log_rerun(self._gui._sess_sel.log.str)


    def _on_right_click(self, event):
        from .gui_table import GUITablePopup
        self.PopupMenu(GUITablePopup(self._gui, self, event,
                                     ['Re-run', 'Save', 'Close'],
                                     ['rerun', 'save', 'close_sess']))

    def _on_close_sess(self, event):
        self._gui._menu_file._on_save(event)
        i = self._gui._sess_list.index(self._gui._sess_sel)
        self._gui._panel_sess._tab.DeleteItem(i)
        del self._gui._sess_list[i]
        del self._gui._sess_item_list[i]
        self._gui._sess_sel = self._gui._sess_list[0]
        self._items = [0]
        self._gui._sess_items = [self._gui._sess_list[0]]
        self._gui._refresh()


    def _on_save(self, event):
        self._gui._menu_file._on_save(event)

    def _on_select(self, event):
        self._sel = event.GetIndex()
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_item_sel.append(self._sel)
        self._entry_select()



    def _on_veto(self, event):
        if event.GetColumn() in [0,2,3,4,5]:
            event.Veto()
        else:
            event.Skip()


    def _refresh(self):
        for i, s in zip(self._items, self._gui._sess_items):
            x = s.spec._safe(s.spec.x)
            self._tab.SetItem(i, 2, "[%3.2f, %3.2f] %s"
                              % (x[0].value, x[-1].value, x.unit))
            self._tab.SetItem(i, 3, str(len(x)))
            try:
                x = s.lines._safe(s.lines.x)
                self._tab.SetItem(i, 4, str(len(x)))
            except:
                pass
            try:
                x = s.systs.z
                self._tab.SetItem(i, 5, str(len(x)))
            except:
                pass

    def _select(self, _sel=0):
        _sel = int(_sel)

        evt = wx.ListEvent()
        evt.SetIndex(_sel)
        self._on_select(evt)


    def json_load(self, path='.'):
        """@brief Load from JSON
        @details Load a set menu from a JSON file.
        @param path Path to file
        @return 0
        """

        logging.info("I'm loading JSON file %s..." % path)

        with open(path) as json_file:
            log = json_file.read()
            log = log.replace('“', '"')
            log = log.replace('”', '"')
            log = log.replace('—', '--')
            load = json.loads(log)

            self._gui._log_run(load)
