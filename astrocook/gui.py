from . import * #version
from .gui_graph import *
from .gui_image import *
from .gui_menu import *
from .gui_table import *
from .message import *
from astropy import table as at
from copy import deepcopy as dc
import json
import logging
import numpy as np
from sphinx.util import docstrings as ds
import wx
import wx.lib.mixins.listctrl as listmix

class GUI(object):
    """ Class for the GUI. """

    def __init__(self, paths=None):
        """ Constructor """

        try:
            print("┌───────────────────┐")
            print("│ ASTROCOOK 🍪  v%3s │" % version)
            print("└───────────────────┘")
        except:
            print("-----------------")
            print(" ASTROCOOK  v%3s " % version)
            print("-----------------")
        print("Cupani et al. 2017-2019 * INAF-OATs")
        self._sess_list = []
        self._sess_sel = None
        self._sess_item_sel = []
        self._menu_spec_id = []
        self._menu_y_conv_id = []
        self._menu_lines_id = []
        self._menu_cont_id = []
        self._menu_nodes_id = []
        self._menu_systs_id = []
        self._menu_mods_id = []
        self._panel_sess = GUIPanelSession(self)
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
                self._panel_sess._on_open(p)

    def _refresh(self, autolim=True, init_cursor=False):
        """ Refresh the GUI after an action """

        self._panel_sess._refresh()
        self._panel_sess._menu._refresh()
        xlim = self._graph_main._graph._ax.get_xlim()
        ylim = self._graph_main._graph._ax.get_ylim()
        goodlim = True
        if xlim == (0.0, 1.0) and ylim == (0.0, 1.0):
            goodlim = False
        if autolim and goodlim and 1==0:
            self._graph_main._refresh(self._sess_items, xlim=xlim)
        else:
            self._graph_main._refresh(self._sess_items)
        if hasattr(self, '_graph_det'):
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
        for s in ['spec', 'lines', 'systs']:
            if hasattr(self, '_tab_'+s):
                if hasattr(getattr(self, '_tab_'+s), '_data'):
                    getattr(self, '_tab_'+s)._on_view(event=None,
                                                      from_scratch=False)
                if hasattr(self, '_col_sel') \
                    and self._col_sel < self._col_tab.GetNumberCols():
                    self._col_values = \
                        [float(self._col_tab.GetCellValue(i, self._col_sel)) \
                         for i in range(self._col_tab.GetNumberRows())]

        if hasattr(self, '_graph_hist'):
            self._graph_hist._refresh(self._sess_items)


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

        # Import GUI
        self._gui = gui
        self._gui._panel_sess = self

        # Create table
        panel = wx.Panel(self)
        self._tab = GUIControlList(panel, 0)
        self._tab.InsertColumn(0, "name", width=200)
        self._tab.InsertColumn(1, "object", width=150)
        self._tab.InsertColumn(2, "active range", width=200)
        self._tab.InsertColumn(3, "# rows", width=100)
        self._tab.InsertColumn(4, "# nodes", width=100)
        self._tab.InsertColumn(5, "# lines", width=100)
        self._tab.InsertColumn(6, "# systems", width=100)
        self._tab.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self._on_veto)
        self._tab.Bind(wx.EVT_LIST_END_LABEL_EDIT, self._on_edit)
        self._tab.Bind(wx.EVT_LIST_ITEM_SELECTED, self._on_select)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._tab, 1, wx.EXPAND)
        panel.SetSizer(self._box)
        self._menu = GUIMenu(self._gui)
        self.SetMenuBar(self._menu.bar())
        self.Show()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _on_add(self, sess, open=True):
        # _sel is the last selection; _items is the list of all selections.
        self._sel = self._tab.GetItemCount()
        self._items = [self._sel]

        self._tab._insert_string_item(self._sel, "%s (%s)"
                                     % (sess.name, str(self._sel)))
        self._gui._sess_list.append(sess)

        # Similarly, _sess_sel contains the last selected session; _sess_items
        # contains all selected sessions
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_items = [self._gui._sess_sel]
        if open:
            self._gui._sess_sel.open()
        x = sess.spec._safe(sess.spec.x)#.value
        self._gui._refresh(autolim=False)

        # Enable import from depending on how many sessions are present
        edit = self._menu._edit
        edit._menu.Enable(edit._start_id+300, len(self._gui._sess_list)>0)
        edit._menu.Enable(edit._start_id+301, len(self._gui._sess_list)>1)


    def _on_edit(self, event):
        self._gui._sess_list[self._sel].spec.meta['object'] = event.GetLabel()


    def _on_open(self, path):
        """ Behaviour for Session > Open """

        #name = path.split('/')[-1][:-5]
        name = path.split('/')[-1].split('.')[0]
        logging.info("I'm loading session %s..." % path)
        sess = Session(path=path, name=name)
        self._gui._panel_sess._on_add(sess, open=True)
        if sess._open_twin:
            logging.info("I'm loading twin session %s..." % path)
            sess = Session(path=path, name=name, twin=True)
            self._gui._panel_sess._on_add(sess, open=True)


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
        exit()

    def _on_select(self, event):
        self._sel = event.GetIndex()
        self._gui._sess_sel = self._gui._sess_list[self._sel]
        self._gui._sess_item_sel = self._tab._get_selected_items()

        # Enable session combine depending on how many sessions are selected
        edit = self._menu._edit
        edit._menu.Enable(edit._start_id+302, len(self._gui._sess_item_sel)>1)

        item = self._tab.GetFirstSelected()
        self._items = []
        while item != -1:
            self._items.append(item)
            item = self._tab.GetNextSelected(item)
        self._gui._sess_items = [self._gui._sess_list[i] for i in self._items]
        self._gui._refresh()

    def _on_veto(self, event):
        if event.GetColumn() in [0,2,3,4,5]:
            event.Veto()
        else:
            event.Skip()

    def _refresh(self):
        for i, s in zip(self._items, self._gui._sess_items):
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
            logging.error("Attribute %s is None." % attrn)
            return None


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


    def combine(self, name='*_combined'):
        """ @brief Combine two or more sessions
        @details When sessions are combined, a new session is created, with a
        new spectrum containing all entries from the spectra of the combined
        sessions. Other objects from the sessions (line lists, etc.) are
        discarded.
        @param name Name of the output session
        @return Combined session
        """
        name_in = name
        #sel = self._tab._get_selected_items()
        sel = self._gui._sess_item_sel
        spec = dc(self._gui._sess_list[sel[0]].spec)
        if name_in[0] == '*':
            name = self._gui._sess_list[sel[0]].name
        for s in sel[1:]:
            spec._t = at.vstack([spec._t, self._gui._sess_list[s].spec._t])
            if name_in[0] == '*':
                name += '_' + self._gui._sess_list[s].name
        if name_in[0] == '*':
            name += name_in[1:]
        sess = Session(name=name, spec=spec)
        return sess

    def load_json(self, path='.'):
        """@brief Load from JSON
        @details Load a set menu from a JSON file.
        @param path Path to file
        @return 0
        """

        with open(path) as json_file:
            d = json.load(json_file)
            for r in d['set_menu']:
                if r['cookbook'][:8]=='cookbook' or r['cookbook']=='cb':
                    cb = self._gui._sess_sel.cb
                elif r['cookbook'] == '':
                    cb = self._gui
                else:
                    rs = r['cookbook'].split('.')
                    cb = getattr(self._gui, rs[0])
                    for s in rs[1:]:
                        cb = getattr(cb, s)
                getattr(cb, r['recipe'])(**r['params'])


    def struct_compare(self, struct_A='0,spec,x', struct_B='0,spec,y',
                       struct_out='0,spec,diff', op='subtract'):


    #struct_A='0,systs,z', struct_B='1,systs,z',
                       #struct_out='0,systs,diff_z', op='subtract'):
        """ @brief Compare structures
        @details Compare two structures from the same session or two different
        sessions.
        @param struct_A Structure A (session,table,column)
        @param struct_B Structure B (same syntax) or scalar
        @param struct_out Output structure (same syntax)
        @param op Binary operator
        @return 0
        """

        parse_A = self._struct_parse(struct_A, length=3)
        parse_out = self._struct_parse(struct_out, length=2)
        if parse_A is None or parse_out is None: return 0
        coln_A, col_A, _ = parse_A
        attrn_out, attr_out, all_out = parse_out
        try:
            col_B = np.full(np.shape(col_A), float(struct_B))
        except:
            parse_B = self._struct_parse(struct_B, length=3)
            parse_out = self._struct_parse(struct_out, length=2)
            if parse_B is None: return 0
            coln_B, col_B, _ = parse_B
            if len(col_A) != len(col_B):
                logging.error("The two columns have different lengths! %s" \
                            % msg_try_again)
                return 0

        if len(col_A) != len(attr_out._t):
            logging.error("The output table have different length than the "
                          "input columns! %s" \
                          % msg_try_again)
            return 0

        if not hasattr(np, op):
            logging.error("Numpy doesn't have a %s operator." % op)
            return 0

        getattr(self._gui._sess_list[all_out[0]], all_out[1])._t[all_out[2]] = \
            getattr(np, op)(col_A, col_B)

        return 0


    def struct_import(self, struct='0,systs', mode='replace'):
        """ @brief Import structure
        @details Import a structure from a session into the current one. The
        structure is either replaced or appended to the corresponding one.
        @param struct Structure (session number,table)
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
            setattr(self._gui._sess_sel, attrn, attr)

        if mode=='append':
            getattr(self._gui._sess_sel, attrn)._append(dc(attr))

        if attrn=='systs':
            self._gui._sess_sel.cb._spec_update()

        return 0
