from . import * #version
from .defaults import *
from .functions import expr_eval
from .gui_graph import *
from .gui_image import *
from .gui_log import *
from .gui_menu import *
from .gui_table import *
from .message import *
from astropy import table as at
from collections import OrderedDict
from copy import deepcopy as dc
import json
import logging
from matplotlib import pyplot as plt
import numpy as np
from sphinx.util import docstrings as ds
import wx
import wx.lib.mixins.listctrl as listmix

class GUI(object):
    """ Class for the GUI. """

    def __init__(self, paths=None):
        """ Constructor """

        try:
            l = ['─']*(16+len(version))
            print("┌%s┐" % ''.join(l))
            print("│ ASTROCOOK 🍪 v%3s │" % version)
            print("└%s┘" % ''.join(l))
        except:
            l = ['-']*(17+len(version))
            print(''.join(l))
            print(" ASTROCOOK  v%3s " % version)
            print(''.join(l))
        print("Cupani et al. 2017-2020 * INAF-OATs")
        self._sess_list = []
        self._sess_item_list = []
        #self._graph_elem_list = []
        #self._meta_list = []
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
        self._menu_tab_id = []
        self._panel_sess = GUIPanelSession(self)
        self._id_zoom = 9
        self._data_lim = None
        self._tag = ""
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

    def _log_rerun(self, log):
        i = self._sess_list.index(self._sess_sel)
        self._panel_sess._tab.DeleteItem(i)
        del self._sess_list[i]
        del self._sess_item_list[i]

        # Run selected log
        self._log_run(json.loads(log), skip_tab=True)


    def _log_run(self, load, skip_tab=False):

#        for obj in ['spec', 'lines', 'systs']:
#            if hasattr(self, '_tab_'+obj) and hasattr(getattr(self, '_tab_'+obj), '_data'):
#                tab = getattr(self, '_tab_'+obj)
#                tab._on_close(None)

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
            #if hasattr(cb, '_menu_view'): print(getattr(self, '_tab_'+r['params']['obj']+'_shown'))
            if hasattr(cb, '_menu_view'):
                try:
                    if getattr(self, '_tab_'+r['params']['obj']+'_shown') \
                        and r['params']['check']==True:
                        skip = True
                except:
                    pass
            #print(r['recipe'], skip)

            """
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

                if r['recipe']=='_data_init' and getattr(self, tab)._shown:
                    skip = True
            """
            if not skip: out = getattr(cb, r['recipe'])(**r['params'])
            if out is not None and out != 0:
                self._panel_sess._on_add(out, open=False)

            #print(r['recipe'])
            if out is None or out==0:
                #print(cb.__dict__)
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


    def _refresh(self, init_cursor=False, init_tab=True, autolim=True,
                 autosort=True, _xlim=None, _ylim=None):
        """ Refresh the GUI after an action """


        self._panel_sess._refresh()
        self._panel_sess._menu._refresh()
        if hasattr(self, '_dlg_mini_defs') \
            and self._dlg_mini_defs._shown:
            self._dlg_mini_defs._refresh()
        if hasattr(self, '_dlg_mini_graph') \
            and self._dlg_mini_graph._shown:
            self._graph_main._elem = self._sess_sel._graph_elem
            self._graph_main._lim = self._sess_sel._graph_lim
            self._dlg_mini_graph._refresh()
        else:
            if hasattr(self._sess_sel, '_graph_elem'):
                self._graph_main._elem = self._sess_sel._graph_elem
            else:
                self._graph_main._elem = elem_expand(graph_elem,
                    self._panel_sess._sel)
            if hasattr(self._sess_sel, '_graph_lim'):
                self._graph_main._lim = self._sess_sel._graph_lim
            else:
                self._graph_main._lim = graph_lim_def
            """
            try:
                if hasattr(self._sess_sel, '_graph_elem'):
                    self._graph_det._elem = self._sess_sel._graph_elem
                else:
                    self._graph_det._elem = elem_expand(graph_elem,
                        self._panel_sess._sel)
            except:
                pass
            """

        if hasattr(self, '_dlg_mini_log') and self._dlg_mini_log._shown:
            self._dlg_mini_log._refresh()
        if hasattr(self, '_dlg_mini_meta') and self._dlg_mini_meta._shown:
            self._dlg_mini_meta._refresh()
        """
        else:
            self._dlg_mini_meta._meta = meta_parse(self._sess_sel.spec._meta)
            self._dlg_mini_meta._meta_backup = meta_parse(
                self._sess_sel.spec._meta_backup)
        """

        ax = self._graph_main._graph._ax
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        dl = self._data_lim
        """
        try:
            print(xlim[0], dl[0], xlim[1], dl[1])
            print(ylim[0], dl[2], ylim[1], dl[3])
        except:
            pass
        """
        if dl is None or ((xlim[0] <= dl[0] or dl[0]==0.) \
                      and (xlim[1] >= dl[1] or dl[1]==1.) \
                      and (ylim[0] <= dl[2] or dl[2]==0.) \
                      and (ylim[1] >= dl[3] or dl[3]==1.)):
            self._graph_main._graph._zoom = False
            dl = (xlim[0], xlim[1], ylim[0], ylim[1])
        #print(self._graph_main._graph._zoom)


        #axc.set_xlim(0, 1)
        #axc.set_ylim(0, 1)
        #print(self._graph_main._graph._zoom)
        #if xlim == axc.get_xlim() and ylim == axc.get_ylim():
        #    print('false')
        #    self._graph_main._graph._zoom = False



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

        if hasattr(self, '_graph_det'):
            #self._refresh_graph_det(init_cursor=init_cursor, autolim=autolim)
            #"""
            #print(self._sess_sel.__dict__)
            if self._sess_sel.systs is None:
                self._graph_det._on_close()
            else:
                graph = self._graph_det._graph
                if hasattr(graph, '_axes'):
                    for key in graph._zems:
                        xunit = self._sess_sel.spec.x.unit
                        #print('before', xunit)
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
                        #print('after', xunit)
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
                #"""
        for s in ['spec', 'lines', 'systs']:
            if hasattr(self, '_tab_'+s) and init_tab:
                if hasattr(getattr(self, '_tab_'+s), '_data'):
                    #print(getattr(self._sess_sel, s))
                    #print(getattr(self, '_tab_'+s))
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
                    self._col_values = \
                        [float(self._col_tab.GetCellValue(i, self._col_sel)) \
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
        """
        self._tab.InsertColumn(0, "name", width=200)
        self._tab.InsertColumn(1, "object", width=150)
        self._tab.InsertColumn(2, "active range", width=200)
        self._tab.InsertColumn(3, "# rows", width=100)
        self._tab.InsertColumn(4, "# nodes", width=100)
        self._tab.InsertColumn(5, "# lines", width=100)
        self._tab.InsertColumn(6, "# systems", width=100)
        """
        self._tab.InsertColumn(0, "name", width=350)
        self._tab.InsertColumn(1, "active range", width=200)
        self._tab.InsertColumn(2, "# rows", width=100)
        self._tab.InsertColumn(3, "# lines", width=100)
        self._tab.InsertColumn(4, "# systems", width=100)
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
        #print(self._gui._sess_item_list)

        sess.log = GUILog(self._gui)
        sess.defs = Defaults(self._gui)

        missing = []
        for i in range(self._tab.GetItemCount()+1):
            if i not in self._gui._sess_item_list:
                missing.append(i)

        self._sel = missing[0]
        #print(self._sel)
        self._items = [self._sel]
        self._tab._insert_string_item(self._sel, "%s (%s)"
                                     % (sess.name, str(self._sel)))
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
            self._gui._sess_sel.open()

        x = sess.spec._safe(sess.spec.x)#.value
        #self._gui._graph_elem_list.append(self._gui._graph_main._elem)
        self._gui._sess_sel._graph_elem = elem_expand(graph_elem, self._sel)
        self._gui._sess_sel._graph_lim = graph_lim_def
        #print(self._gui._sess_sel._graph_elem)
        #self._gui._meta_list.append(self._gui._dlg_mini_meta._meta)
        #self._gui._refresh(autolim=False)

        self._gui._refresh(init_tab=False, autolim=False)

        # Enable import from depending on how many sessions are present
        edit = self._menu._edit
        #edit._menu.Enable(edit._start_id+300, len(self._gui._sess_list)==2)
        #edit._menu.Enable(edit._start_id+301, len(self._gui._sess_list)>1)
        edit._menu.Enable(edit._start_id+310, len(self._gui._sess_list)>0)
        edit._menu.Enable(edit._start_id+311, len(self._gui._sess_list)>0)



    def _on_edit(self, event):
        self._gui._sess_list[self._sel].spec.meta['object'] = event.GetLabel()


    def _on_open(self, path):
        """ Behaviour for Session > Open """

        name = '.'.join(path.split('/')[-1].split('.')[:-1])
        logging.info("I'm loading session %s..." % path)
        sess = Session(gui=self._gui, path=path, name=name)
        self._gui._panel_sess._on_add(sess, open=True)

        if sess._open_twin:
            logging.info("I'm loading twin session %s..." % path)
            sess = Session(gui=self._gui, path=path, name=name, twin=True)
            self._gui._panel_sess._on_add(sess, open=True)

        sess.log.append_full('_panel_sess', '_on_open', {'path': path})


    def _entry_select(self):
        try:
            self._gui._sess_sel.cb.sess = self._gui._sess_sel
        except:
            pass

        # Enable session equalize/combine depending on how many sessions are selected
        edit = self._menu._edit
        edit._menu.Enable(edit._start_id+300, len(self._gui._sess_item_sel)==2)
        edit._menu.Enable(edit._start_id+301, len(self._gui._sess_item_sel)>1)

        item = self._tab.GetFirstSelected()
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
        """
        self._gui._panel_sess.Close()
        self._gui._graph_main.Close()
        self._gui._tab_spec.Close()
        self._gui._tab_lines.Close()
        self._gui._tab_systs.Close()
        self._gui._tab_mods.Close()
        """
        #exit()
        os._exit(1)

    def _on_deselect(self, event):
        self._sel = event.GetIndex()
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_item_sel = []
        self._entry_select()


    def _on_rerun(self, event):
        from .gui_dialog import GUIDialogMiniLog
        self._gui._log_rerun(self._gui._sess_sel.log.str)


    def _on_right_click(self, event):
        from .gui_table import GUITablePopup
        self.PopupMenu(GUITablePopup(self._gui, self, event,
                                     ['Re-run', 'Save'],
                                     ['rerun', 'save']))


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
            """
            obj = s.spec.meta['object']
            self._tab.SetItem(i, 1, obj)
            x = s.spec._safe(s.spec.x)
            self._tab.SetItem(i, 2, "[%3.2f, %3.2f] %s"
                              % (x[0].value, x[-1].value, x.unit))
            self._tab.SetItem(i, 3, str(len(x)))
            try:
                x = s.nodes._safe(s.nodes.x)
                self._tab.SetItem(i, 4, str(len(x)))
            except:
                pass
            try:
                x = s.lines._safe(s.lines.x)
                self._tab.SetItem(i, 5, str(len(x)))
            except:
                pass
            try:
                x = s.systs.z
                self._tab.SetItem(i, 6, str(len(x)))
            except:
                pass
            """
            x = s.spec._safe(s.spec.x)
            self._tab.SetItem(i, 1, "[%3.2f, %3.2f] %s"
                              % (x[0].value, x[-1].value, x.unit))
            self._tab.SetItem(i, 2, str(len(x)))
            try:
                x = s.lines._safe(s.lines.x)
                self._tab.SetItem(i, 3, str(len(x)))
            except:
                pass
            try:
                x = s.systs.z
                self._tab.SetItem(i, 4, str(len(x)))
            except:
                pass

    def _select(self, _sel=0):

        _sel = int(_sel)
        """
        sel = self._gui._sess_item_sel
        sess_list = self._gui._sess_list
        if sel == []:
            try:
                sel = [int(s) \
                       for s in _sel.replace('[','').replace(']','').split(',')]
            except:
                pass
        if sel == []:
            sel = range(len(sess_list))
        """
        evt = wx.ListEvent()
        evt.SetIndex(_sel)
        self._on_select(evt)


    def _struct_parse(self, struct, length=2):
        sess_list = self._gui._sess_list

        parse = struct.split(',')

        if len(parse) < length:
            logging.error("I can't parse the structure.")
            return None

        # Session
        sessn = parse[0]
        try:
            sessn = int(sessn)
            parse[0] = sessn
        except:
            logging.error(msg_param_fail)
            return None

        if sessn > len(sess_list):
            logging.error("I can't find session %s." % sessn)
            return None

        # Attribute
        attrn = parse[1]
        sess = sess_list[sessn]
        if not hasattr(sess, attrn):
            logging.error(msg_attr_miss(attrn))
            return None

        attr = getattr(sess, attrn)
        if attr is None:
            logging.warning("Attribute %s is None." % attrn)
            return attrn, attr, parse


        if length==3:
            # Column
            coln = parse[2]
            if coln not in getattr(sess, attrn)._t.colnames:
                logging.error(msg_col_miss(coln))
                return None
            col = getattr(sess, attrn)._t[coln]
            return coln, col, parse
        else:
            return attrn, attr, parse


    def combine(self, name='*_combined', _sel=''):
        """ @brief Combine two or more sessions
        @details Combine two or more sessions. A new session is created, with a
        new spectrum containing all entries from the spectra of the combined
        sessions. Other objects from the sessions (line lists, etc.) are
        discarded.
        @param name Name of the output session
        @return Combined session
        """
        name_in = name
        #sel = self._tab._get_selected_items()
        sel = self._gui._sess_item_sel
        sess_list = self._gui._sess_list

        """
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            try:
                sel = [int(s) \
                       for s in _sel.replace('[','').replace(']','').split(',')]
            except:
                pass
        if sel == []:
            sel = range(len(sess_list))
        self._gui._sess_item_sel = sel
        """
        sel = _sel

        struct_out = {}
        for struct in sess_list[sel[0]].seq:
            struct_out[struct] = dc(getattr(sess_list[sel[0]], struct))


        if name_in[0] == '*':
            name = sess_list[sel[0]].name

        logging.info("Combining sessions %s..." % ', '.join(str(s) for s in sel))
        for s in sel[1:]:
            #spec._t = at.vstack([spec._t, self._gui._sess_list[s].spec._t])

            for struct in sess_list[s].seq:
                if getattr(sess_list[s], struct) != None:
                    if struct_out[struct] != None:
                        struct_out[struct]._append(
                            getattr(sess_list[s], struct))
                    else:
                        struct_out[struct] = dc(getattr(sess_list[s], struct))

            if name_in[0] == '*':
                name += '_' + sess_list[s].name

        struct_out['spec']._t.sort('x')
        if name_in[0] == '*':
            name += name_in[1:]
        sess = Session(gui=self._gui, name=name, spec=struct_out['spec'],
                       nodes=struct_out['nodes'], lines=struct_out['lines'],
                       systs=struct_out['systs'])
        return sess


    def equalize(self, xmin, xmax, _sel=''):
        """ @brief Equalize two sessions
        @details Equalize the flux level of one session to another one. The
        last-selected session is equalized to the first-selected one. The
        equalization factor is the ratio of the median flux within the
        specified wavelength interval.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return 0
        """

        try:
            xmin = float(xmin) * au.nm
            xmax = float(xmax) * au.nm
        except ValueError:
            logging.error(msg_param_fail)
            return None

        """
        sel = self._gui._sess_item_sel
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            sel = [int(s) \
                for s in _sel.replace('[','').replace(']','').split(',')]
        self._gui._sess_item_sel = sel
        """
        sel = _sel
        logging.info("Equalizing session %i to session %i... "
                     % (sel[1], sel[0]))

        for i,s in enumerate(sel):
            sess = self._gui._sess_list[s]
            w = np.where(np.logical_and(sess.spec.x>xmin, sess.spec.x<xmax))[0]
            if len(w)==0:
                logging.error("I can't use this wavelength range for "
                              "equalization. Please choose a range covered by "
                              "both sessions.")
                return(0)
            if i == 0:
                f = np.nanmedian(sess.spec.y[w]).value
                #print(np.median(sess.spec.y[w]))
            else:
                f = f/np.nanmedian(sess.spec.y[w]).value
                #print(np.median(sess.spec.y[w]), f)
                logging.info("Equalization factor: %3.4f." % f)
                sess.spec.y = f*sess.spec.y
                sess.spec.dy = f*sess.spec.dy

        return 0


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


    def struct_modify(self, col_A='0,spec,x', col_B='0,spec,y',
                      col_out='0,spec,diff', op='subtract'):
        """ @brief Modify a data structure using a binary operator
        @details Modify a data structure using a binary operator. An output
        column is computed applying a binary operator to two input columns, or
        an input column and a scalar. Columns are described by a string with the
        session number, the structure tag (spec, lines, systs), and the column
        name separated by a comma (e.g. 0,spec,x, meaning "column x of spectrum
        from session 0"). They can be from different data structures only if
        they have the same length. If the output column already exists, it is
        overwritten.
        @param col_A Structure A
        @param col_B Structure B or scalar
        @param col_out Output structure
        @param op Binary operator
        @return 0
        """

        parse_A = self._struct_parse(col_A, length=3)
        parse_out = self._struct_parse(col_out, length=2)
        #print(parse_A, parse_out)
        if parse_A is None or parse_out is None: return 0
        coln_A, colp_A, _ = parse_A
        attrn_out, attr_out, all_out = parse_out
        try:
            colp_B = np.full(np.shape(colp_A), float(col_B))
        except:
            parse_B = self._struct_parse(col_B, length=3)
            parse_out = self._struct_parse(col_out, length=2)
            if parse_B is None: return 0
            coln_B, colp_B, _ = parse_B
            if len(colp_A) != len(colp_B):
                logging.error("The two columns have different lengths! %s" \
                            % msg_try_again)
                return 0

        if len(colp_A) != len(attr_out._t):
            logging.error("The output table have different length than the "
                          "input columns! %s" \
                          % msg_try_again)
            return 0

        if not hasattr(np, op):
            logging.error("Numpy doesn't have a %s operator." % op)
            return 0
        getattr(self._gui._sess_list[all_out[0]], all_out[1])._t[all_out[2]] = \
            getattr(np, op)(colp_A, colp_B)

        return 0


    def struct_modify2(self, col='', expr=''):
        """ @brief Modify a data structure using a binary operator
        @details Modify a data structure using a binary operator. An output
        column is computed from an expression with input columns as arguments.
        The expression must be parsable by AST, with columns described by a
        string with the session number, the structure tag (spec, lines, systs),
        and the column name separated by a comma (e.g. 0,spec,x, meaning "column
        x of spectrum from session 0"). Columns can be from different data
        structures only if they have the same length. If the output column
        already exists, it is overwritten.
        @param col Output column
        @param expr Expression
        @return 0
        """

        # Reversed to parse sessions with higher number first, and avoid overwriting
        sess_list = self._gui._sess_list[::-1]
        for i, s in enumerate(sess_list):
            if s.spec is not None:
                for c in sorted(s.spec._t.colnames, key=len, reverse=True):
                    expr = expr.replace('%i,spec,%s' % (len(sess_list)-1-i, c),
                                        str(list(np.array(s.spec._t[c]))))
            if s.lines is not None:
                for c in sorted(s.lines._t.colnames, key=len, reverse=True):
                    expr = expr.replace('%i,spec,%s' % (len(sess_list)-1-i, c),
                                        str(list(np.array(s.lines._t[c]))))
            if s.systs is not None:
                for c in sorted(s.systs._t.colnames, key=len, reverse=True):
                    expr = expr.replace('%i,spec,%s' % (len(sess_list)-1-i, c),
                                        str(list(np.array(s.systs._t[c]))))

        #print(expr)
        #print(len(expr))
        out = expr_eval(ast.parse(expr, mode='eval').body)

        _, _, all_out = self._struct_parse(col, length=2)
        struct = getattr(self._gui._sess_list[all_out[0]], all_out[1])
        if all_out[2] in struct._t.colnames: # and False:
            col_out = struct._t[all_out[2]]
            try:
                struct._t[all_out[2]] = expr_eval(ast.parse(expr, mode='eval').body) * col_out.unit
            except:
                struct._t[all_out[2]] = expr_eval(ast.parse(expr, mode='eval').body)
        else:
            struct._t[all_out[2]] = expr_eval(ast.parse(expr, mode='eval').body)

        return 0


    def struct_import(self, struct='0,systs', mode='replace'):
        """ @brief Import a data structure from a session into the current one
        @details The structure to be imported is described by a string with the
        session number and the structure tag (spec, lines, systs), separated by
        a comma (e.g. 0,spec, meaning "spectrum from session 0"). The imported
        structure is either replaced or appended to the corresponding one in the
        current session.
        @param struct Structure
        @param mode Mode (replace or append)
        @return 0
        """

        parse = self._struct_parse(struct)
        if parse is None: return 0
        attrn, attr, _ = parse
        attr = dc(attr)

        if attrn == 'systs' \
            and 'cont' not in self._gui._sess_sel.spec.t.colnames:
            logging.error("Attribute %s requires a continuum. Please try "
                          "Recipes > Guess continuum before." % attrn)
            return 0

        if mode=='replace':
            if attr is None:
                logging.warning("I'm replacing structure with None.")
                setattr(self._gui._sess_sel, attrn, attr)
                return 0
            if attrn in ['lines', 'systs']:
                #spec = self._gui._sess_sel.spec
                x = self._gui._sess_sel.spec.x.to(au.nm)
                attr = attr._region_extract(np.min(x), np.max(x))

                # Redefine regions from spectrum
            if attrn == 'systs':
                for m in attr._mods_t:
                    mod = m['mod']
                    mod._spec = self._gui._sess_sel.spec
                    #mod._xf, mod._yf, mod._wf, mod._ys = \
                    mod._xf, mod._yf, mod._wf = \
                        mod._make_regions(mod, mod._spec._safe(mod._spec.x)\
                                               .to(au.nm).value)
            setattr(self._gui._sess_sel, attrn, attr)

        if mode=='append':
            if attr is None:
                logging.warning("I'm not appending None.")
                return 0
            attr_dc = dc(attr)
            if attrn == 'systs':
                id_max = np.max(getattr(self._gui._sess_sel, attrn)._t['id'])
                attr_dc._t['id'] = attr_dc._t['id']+id_max
            #print(len(attr_dc._t))
            #print(len(np.unique(attr_dc._t['id'])))
            #print(len(attr_dc._t))
            #print(len(getattr(self._gui._sess_sel, attrn)._t))
            getattr(self._gui._sess_sel, attrn)._append(attr_dc)

        if attrn=='systs':
            self._gui._sess_sel.cb._mods_recreate()
            self._gui._sess_sel.cb._spec_update()

        return 0
