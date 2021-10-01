import astropy.units as au
from .functions import elem_expand, meta_parse, trans_parse
from .message import *
from .vars import graph_elem, graph_lim_def, hwin_def, xem_d
from collections import OrderedDict
from copy import deepcopy as dc
import inspect
import json
import numpy as np
import datetime as dt
import wx

class GUIDialog(wx.Dialog):

    def __init__(self,
                 gui,
                 title,
                 attr,
                 obj=None):
        self._gui = gui
        self._gui._dlg = self
        super(GUIDialog, self).__init__(parent=None, title=title)
        self._attr = np.array(attr, ndmin=1)
        self._obj = obj
        self._methods = []
        self._params = []
        self._brief = []
        self._details = []
        self._doc = []
        self._ctrl = []

        for i, a in enumerate(self._attr):
            self._sess_sel =  self._gui._sess_sel
            self._cb =  self._gui._sess_sel.cb
            if self._obj == None:
                if hasattr(self._cb, a):
                    self._obj = self._cb
                else:
                    self._obj = self._sess_sel
            method = getattr(self._obj, a)
            self._methods.append(method)
            self._get_params(method)
            self._get_last(method)
            self._get_doc(method)

        self._panel = wx.Panel(self)
        self._bottom = wx.BoxSizer(wx.VERTICAL)
        self._core = wx.BoxSizer(wx.VERTICAL)

    def _box_buttons(self, cancel_run=True):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        if cancel_run:
            cancel_button = wx.Button(self, label='Cancel')
            cancel_button.Bind(wx.EVT_BUTTON, self._on_cancel)
            run_button = wx.Button(self, label='Run')
            run_button.Bind(wx.EVT_BUTTON, self._on_run)
            run_button.SetDefault()
            buttons.Add(cancel_button, 0, wx.RIGHT, border=5)
            buttons.Add(run_button, 0)
        else:
            ok_button = wx.Button(self, label='OK')
            ok_button.Bind(wx.EVT_BUTTON, self._on_ok)
            ok_button.SetDefault()
            buttons.Add(ok_button, 0, border=5)

        self._bottom.Add(self._panel, 0, wx.EXPAND|wx.ALL, border=10)
        self._bottom.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        self._bottom.SetSizeHints(self)

    def _get_doc(self, method):
        full = inspect.getdoc(method)
        split = full.split('@')
        self._brief.append([s[6:-1] for s in split \
                           if s[0:5]=='brief'][0].replace('\n', ' '))
        self._details.append([s[8:-1] for s in split \
                              if s[0:7]=='details'][0].replace('\n', ' '))
        self._doc.append([s[6:-1].split(' ', 1)[1] \
                          for s in split if s[0:5]=='param'])

    def _get_last(self, method):
        if self._params_last is not None:
            for pls in self._params_last:
                for pl in pls:
                    for p in self._params[0]:
                        if p == pl:
                            self._params[0][p] = pls[pl]


    def _get_params(self, method):
        keys = inspect.getargspec(method)[0][1:]
        defs = inspect.getargspec(method)[-1]
        if defs == None:
            defs = []
        defs = [str(d) for d in defs]
        values = np.append(['']*(len(keys)-len(defs)), defs)
        self._params.append(OrderedDict(zip(keys, values)))

    def _on_cancel(self, e):
        if hasattr(self._gui, '_dlg_mini_systems'):
            self._gui._dlg_mini_systems._cursor_refresh()
            """
            # Unfreeze cursors in case they were frozen
            self._gui._graph_main._graph._cursor_frozen = False
            self._gui._graph_det._graph._cursor_frozen = False

            # Refresh cursor
            shown = True if self._gui._dlg_mini_systems._shown else False
            self._gui._dlg_mini_systems._shown = False
            self._gui._dlg_mini_systems._on_cancel(e=None)
            self._gui._dlg_mini_systems._shown = True
            self._gui._dlg_mini_systems._on_apply(e=None)
            self._gui._dlg_mini_systems._shown = shown
            """
        self.Close()

    def _on_ok(self, e):
        self._update_params()
        for pps in self._params_parent:
            for pp in pps:
                for p in self._params[0]:
                    if p == pp:
                        pps[pp] = self._params[0][p]
        self.Close()

    def _on_run(self, e):
        self._update_params()
        for a, p_l in zip(self._attr, self._params):

            m = getattr(self._obj, a)
            logging.info("I'm launching %s..." % a)

            if '_sel' in p_l:
                p_l['_sel'] = self._gui._sess_item_sel
            """
            try:
                p_l['_sel'] = self._gui._sess_item_sel
            except:
                logging.info('No selected session.')
            """

            start = dt.datetime.now()
            out = m(**p_l)
            end = dt.datetime.now()

            logging.info("I completed %s in %3.3f seconds!" \
                         % (a, (end-start).total_seconds()))

            sel_old = self._gui._sess_list.index(self._gui._sess_sel)

            if out is not None:
                if out is 0:
                    #self._gui._refresh()
                    pass
                else:
                    self._gui._panel_sess._on_add(out, open=False)
                self.Close()

            if out is None or out==0:
                if '_sel' in p_l:
                    sess_list = [self._gui._sess_list[s] for s in p_l['_sel']]
                    self._gui._sess_sel.log.merge_full(self._obj._tag, a, p_l,
                                                       sess_list,
                                                       self._gui._sess_sel)
                else:
                    self._gui._sess_sel.log.append_full(self._obj._tag, a, p_l)
            else:
                if '_sel' in p_l:
                    sess_list = [self._gui._sess_list[s] for s in p_l['_sel']]
                    self._gui._sess_sel.log.merge_full(self._obj._tag, a, p_l,
                                                       sess_list,
                                                       self._gui._sess_sel)
                else:
                    sess_list = [self._gui._sess_list[sel_old]]
                    self._gui._sess_sel.log.merge_full(self._obj._tag, a, p_l,
                                                       sess_list,
                                                       self._gui._sess_sel)

            self._gui._refresh()
            if hasattr(self._gui, '_dlg_mini_systems'):
                self._gui._dlg_mini_systems._cursor_refresh()
                """
                # Unfreeze cursors in case they were frozen
                self._gui._graph_main._graph._cursor_frozen = False
                self._gui._graph_det._graph._cursor_frozen = False

                # Refresh cursor
                shown = True if self._gui._dlg_mini_systems._shown else False
                self._gui._dlg_mini_systems._shown = False
                self._gui._dlg_mini_systems._on_cancel(e=None)
                self._gui._dlg_mini_systems._shown = True
                self._gui._dlg_mini_systems._on_apply(e=None)
                self._gui._dlg_mini_systems._shown = shown
                """


    def _update_params(self):
        for p_l, c_l in zip(self._params, self._ctrl):
            for p, c in zip(p_l, c_l):
                pmod = c.GetValue()
                p_l[p] = pmod


class GUIDialogMethod(GUIDialog):

    def __init__(self,
                 gui,
                 title,
                 attr,
                 obj=None,
                 cancel_run=True,
                 params_last=None,
                 params_parent=None):

        self._params_last = params_last
        self._cancel_run = cancel_run
        self._params_parent = params_parent
        super(GUIDialogMethod, self).__init__(gui, title, attr, obj)
        self._box_descr()
        self._box_params()
        self._box_buttons(self._cancel_run)
        self.SetSizer(self._bottom)
        self.Centre()
        self.SetPosition((self.GetPosition()[0], wx.DisplaySize()[1]*0.25))
        self.Show()

    def _box_descr(self):

        # Description
        #sb = wx.StaticBox(self._panel, label="Description")
        #sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        sbs = wx.StaticBoxSizer(wx.VERTICAL, self._panel, label="Description")

        st = wx.StaticText(sbs.GetStaticBox(), 1, label='')
        st.SetLabel('\n'.join(self._details))
        st.Wrap(400)
        sbs.Add(st, flag=wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND, border=8)
        #st.SetLabel(''.join(self._details))
        self._core.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        self._panel.SetSizer(self._core)


    def _box_params(self):
        # Parameters
        sb = wx.StaticBox(self._panel, label="Parameters")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        len_params = np.sum([len(i) for i in self._params])
        fgs = wx.FlexGridSizer(len_params, 2, 4, 15)
        fgs_add = []
        for p_l, d_l in zip(self._params, self._doc):
            ctrl_l = []
            for p, d in zip(p_l, d_l):
                stat = wx.StaticText(self._panel, -1, label=d+':')
                ctrl = wx.TextCtrl(self._panel, -1, value=str(p_l[p]))
                fgs_add.append((stat, 1, wx.EXPAND))
                fgs_add.append((ctrl, 1, wx.EXPAND))
                ctrl_l.append(ctrl)
            self._ctrl.append(ctrl_l)
        if np.size(self._ctrl) > 0:
            fgs.AddMany(fgs_add)
            sbs.Add(fgs, flag=wx.ALL|wx.EXPAND, border=8)
            self._core.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        self._core.SetMinSize(width=450, height=100)
        self._panel.SetSizer(self._core)


class GUIDialogMethods(GUIDialog):

    def __init__(self,
                 gui,
                 title,
                 attr,
                 obj=None,
                 params_last=None):

        self._params_last = params_last
        super(GUIDialogMethods, self).__init__(gui, title, attr, obj)
        self._box_methods()#panel, core)
        self._box_buttons()
        self.SetSizer(self._bottom)
        self.Centre()
        self.SetPosition((self.GetPosition()[0], wx.DisplaySize()[1]*0.25))
        self.Show()

    def _box_methods(self):
        sb = wx.StaticBox(self._panel, label="Methods")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        bbs = wx.BoxSizer(wx.VERTICAL)
        for a_l, b_l in zip(self._attr, self._brief):
            button = wx.Button(self._panel, label=b_l)
            button.Bind(wx.EVT_BUTTON, lambda e, a=[a_l], b=b_l: \
                        self._on_method(e, a, b))
            bbs.Add(button, 0, wx.BOTTOM|wx.EXPAND, border=5)
        sbs.Add(bbs, flag=wx.LEFT|wx.RIGHT|wx.TOP, border=5)
        self._core.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        self._panel.SetSizer(self._core)

    def _on_method(self, e, attr, brief):
        GUIDialogMethod(self._gui, brief, attr, cancel_run=False,
                        params_last=self._params_last, params_parent=self._params)

class GUIDialogMini(wx.Dialog):
    def __init__(self,
                 gui,
                 title):
        self._gui = gui
        #self._gui._dlg_mini = self
        #self._targ = targ
        #self._series = series
        #self._z = z
        super(GUIDialogMini, self).__init__(parent=None, title=title)
        self._init()

    def _init(self):
        self._panel = wx.Panel(self)
        self._bottom = wx.BoxSizer(wx.VERTICAL)
        self._core = wx.BoxSizer(wx.VERTICAL)
        self._box_ctrl()
        self._box_buttons()
        self.SetSizer(self._bottom)
        self.Centre()
        self.SetPosition((self.GetPosition()[0], wx.DisplaySize()[1]*0.25))
        self.Show()


class GUIDialogMiniDefaults(GUIDialogMini):
    def __init__(self,
                 gui,
                 title):
        self._gui = gui
        self._gui._dlg_mini_defs = self
        self._sel = dc(self._gui._panel_sess._sel)
        self._defs_str = self._gui._sess_sel.defs.str
        self._defs_dict = self._gui._sess_sel.defs.dict
        super(GUIDialogMiniDefaults, self).__init__(gui, title)
        self.Bind(wx.EVT_CLOSE, self._on_cancel)
        self._shown = True

    def _box_ctrl(self):
        fgs = wx.FlexGridSizer(2, 1, 4, 15)
        descr = wx.StaticText(
                    self._panel, -1,
                    label="When clicking on “Apply”, the GUI is refreshed and the\n"
                          "system models are recreated, but not re-fitted (so the\n"
                          "fitting parameters in the system table may not reflect\n"
                          "the changes). Try “fit systems” to re-fit them.")
        self._ctrl_defs = wx.TextCtrl(self._panel, -1, value=self._defs_str,
                                      size=(400, 300), style = wx.TE_MULTILINE)
        #self._ctrl_z = wx.TextCtrl(self._panel, -1, value="%3.7f" % 10, size=(150, -1))
        fgs.AddMany([(self._ctrl_defs, 1, wx.EXPAND), (descr, 1, wx.EXPAND)])
        self._core.Add(fgs, flag=wx.ALL|wx.EXPAND)
        self._panel.SetSizer(self._core)


    def _box_buttons(self):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        apply_button = wx.Button(self, label='Apply')
        apply_button.Bind(wx.EVT_BUTTON, self._on_apply)
        #apply_button.SetDefault()
        buttons.Add(apply_button, 0, wx.RIGHT, border=5)
        default_button = wx.Button(self, label='Load from file')
        default_button.Bind(wx.EVT_BUTTON, self._on_load)
        buttons.Add(default_button)
        self._bottom.Add(self._panel, 0, wx.EXPAND|wx.ALL, border=10)
        self._bottom.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        self._bottom.SetSizeHints(self)

    def _on_apply(self, e=None, refresh=True, log=True):
        """
        defs = self._ctrl_defs.GetValue()
        defs = defs.replace('“', '"')
        defs = defs.replace('”', '"')
        defs = defs.replace('—', '--')
        self._ctrl_defs.SetValue(defs)
        """
        sess = self._gui._sess_sel
        sd = self._gui._sess_sel.defs
        defs_dict = dict(sess.defs.dict)
        self._set(self._ctrl_defs.GetValue())
        for i in sess.defs.dict:
            for k, v in set(sess.defs.dict[i].items()) - set(defs_dict[i].items()):
                for e in sess.defs._extend:
                    if k not in sess.defs._extend[e]:
                        logging.info("I changed parameter %s %s from %s to %s."
                                     % (i, k, str(defs_dict[i][k]),
                                        str(sess.defs.dict[i][k])))
            #diff = {k: sd.dict[i][k] for k, _ \
            #        in set(sd.dict[i].items()) - set(defs_dict[i].items()) }
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_defs', '_on_apply',
                                 {'e': None, 'refresh': refresh})
        if refresh:
            if hasattr(sess, 'systs') and sess.systs is not None:
                sess.cb._mods_recreate2()
            self._gui._refresh(init_cursor=True, init_tab=False)


    def _on_load(self, e=None, path=None, log=True):
        sess = self._gui._sess_sel

        if path is None:
            wildcard = "JSON files (*.json)|*.json"
            with wx.FileDialog(self, "Open file", '.',
                               wildcard=wildcard,
                               style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) \
                               as fileDialog:
                if fileDialog.ShowModal() == wx.ID_CANCEL:
                    return
                path = fileDialog.GetPath()
        sess.defs.open(path)

        if log:
            sess.log.append_full('_dlg_mini_defs', '_on_load',
                                 {'e': None, 'path': path})

        self._refresh()
        if hasattr(sess, 'systs'):
            sess.cb._mods_recreate2()
        self._gui._refresh(init_cursor=True, init_tab=False)


    def _on_cancel(self, e=None, refresh=True, log=True):
        self._shown = False
        self.Destroy()

    def _refresh(self):
        self._ctrl_defs.SetValue(self._defs_str)


    def _set(self, value, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_defs', '_set',
                                 {'value': value, 'log': False})

        defs = value
        defs = defs.replace('“', '"')
        defs = defs.replace('”', '"')
        defs = defs.replace('—', '--')
        self._ctrl_defs.SetValue(defs)
        self._gui._sess_sel.defs.update(defs)
        self._defs_str = self._gui._sess_sel.defs.str
        self._defs_dict = self._gui._sess_sel.defs.dict



class GUIDialogMiniGraph(GUIDialogMini):
    def __init__(self,
                 gui,
                 title,
                 elem=graph_elem):
        self._gui = gui
        self._gui._dlg_mini_graph = self
        self._sel = dc(self._gui._panel_sess._sel)
        #self._elem = '\n'.join([self._sel+','+r for r in elem.split('\n')])
        if hasattr(self._gui._sess_sel, '_graph_elem'):
            self._elem = self._gui._sess_sel._graph_elem
        else:
            self._elem = elem_expand(elem, self._sel)
        self._lim = self._gui._sess_sel._graph_lim
        super(GUIDialogMiniGraph, self).__init__(gui, title)
        self.Bind(wx.EVT_CLOSE, self._on_cancel)
        self._shown = True


    def _box_ctrl(self):
        fgs = wx.FlexGridSizer(4, 1, 4, 15)
        descr_elem = wx.StaticText(
                    self._panel, -1,
                    label="Each line define a graph element as a set of comma-separated\n"
                          "values. Values are: session, table, x column, y column, mask\n"
                          "column (if any), type of graph (plot, step, scatter), line style\n"
                          "or marker symbol, line width or marker size, color, opacity.")
        descr_lim = wx.StaticText(
                    self._panel, -1,
                    label="Adjust limits for the x and y axis of the main graph. Limits\n"
                          "must be two comma-separated values between parentheses.\n"
                          "Values have the same units as in the graph.")
        self._ctrl_elem = wx.TextCtrl(self._panel, -1, value=self._elem,
                                      size=(360, 200), style = wx.TE_MULTILINE)
        self._ctrl_lim = wx.TextCtrl(self._panel, -1, value=self._lim,
                                      size=(360, 38), style = wx.TE_MULTILINE)
        #self._ctrl_z = wx.TextCtrl(self._panel, -1, value="%3.7f" % 10, size=(150, -1))
        fgs.AddMany([(self._ctrl_elem, 1, wx.EXPAND),
                     (descr_elem, 1, wx.EXPAND),
                     (self._ctrl_lim, 1, wx.EXPAND),
                     (descr_lim, 1, wx.EXPAND),
                                  ])
        self._core.Add(fgs, flag=wx.ALL|wx.EXPAND)
        self._panel.SetSizer(self._core)


    def _box_buttons(self):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        apply_button = wx.Button(self, label='Apply')
        apply_button.Bind(wx.EVT_BUTTON, self._on_apply)
        #apply_button.SetDefault()
        buttons.Add(apply_button, 0, wx.RIGHT, border=5)
        default_button = wx.Button(self, label="Back to default")
        default_button.Bind(wx.EVT_BUTTON, self._on_default)
        buttons.Add(default_button)
        self._bottom.Add(self._panel, 0, wx.EXPAND|wx.ALL, border=10)
        self._bottom.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        self._bottom.SetSizeHints(self)


    def _on_apply(self, e=None, refresh=True, log=True):
        # Elements
        self._elem = self._ctrl_elem.GetValue()
        self._set_elem(self._elem)
        self._gui._graph_main._elem = self._elem
        #self._gui._graph_elem_list[self._sel] = self._elem
        self._gui._sess_sel._graph_elem = self._elem
        if hasattr(self._gui, '_graph_det'):
            self._gui._graph_det._elem = elem_expand(graph_elem, self._sel)

        # Limits
        self._lim = self._ctrl_lim.GetValue()
        self._set_lim(self._lim)
        self._gui._graph_main._lim = self._lim
        self._gui._sess_sel._graph_lim = self._lim

        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_graph', '_on_apply',
                                 {'e': None, 'refresh': refresh})
        if refresh:
            self._gui._refresh(init_cursor=True, init_tab=False)


    def _on_default(self, e=None, refresh=True, log=True):
        # Elements
        self._elem = elem_expand(graph_elem, self._sel)
        self._gui._graph_main._elem = self._elem
        #self._gui._graph_elem_list[self._sel] = self._elem
        self._gui._sess_sel._graph_elem = self._elem

        # Limits
        self._lim = graph_lim_def
        self._gui._graph_main._lim= self._lim
        #self._gui._graph_elem_list[self._sel] = self._elem
        self._gui._sess_sel._graph_lim = self._lim

        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_graph', '_on_default',
                                 {'e': None, 'refresh': refresh})
        if refresh: self._gui._refresh(init_cursor=True, init_tab=False, autolim=False)


    def _on_cancel(self, e=None):
        self._shown = False
        self.Destroy()


    def _refresh(self):
        if self._sel != self._gui._panel_sess._sel:
            self._sel = self._gui._panel_sess._sel
            #self._elem = elem_expand(graph_elem, self._sel)
            #self._elem = self._gui._graph_elem_list[self._sel]
        self._elem = self._gui._sess_sel._graph_elem
        self._lim = self._gui._sess_sel._graph_lim

        self._ctrl_elem.SetValue(self._elem)
        self._ctrl_lim.SetValue(self._lim)
        #self._on_apply(refresh=False)


    def _set_elem(self, value, log=True):
        if log:
            sess = self._gui._sess_sel
            i = self._gui._sess_list.index(self._gui._sess_sel)
            value_log = '\nSESS_SEL,'.join(('\n'+value).split('\n%i,' % i))[1:]
            sess.log.append_full('_dlg_mini_graph', '_set_elem',
                                 {'value': value_log, 'log': False})

        self._elem = value
        self._ctrl_elem.SetValue(value)


    def _set_lim(self, value, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_graph', '_set_lim',
                                 {'value': value, 'log': False})

        self._lim = value
        self._ctrl_lim.SetValue(value)


class GUIDialogMiniLog(GUIDialogMini):
    def __init__(self,
                 gui,
                 title):
        self._gui = gui
        self._title = title
        self._gui._dlg_mini_log = self
        self._sel = dc(self._gui._panel_sess._sel)

        self._log = self._gui._sess_sel.log.str
        super(GUIDialogMiniLog, self).__init__(gui, title)
        self.Bind(wx.EVT_CLOSE, self._on_cancel)
        self._shown = False


    def _box_ctrl(self):
        fgs = wx.FlexGridSizer(2, 1, 4, 15)
        #print(self._log)
        self._ctrl_log = wx.TextCtrl(self._panel, -1, value=self._log,
                                     size=(400, 300),
                                     style = wx.TE_MULTILINE)#|wx.TE_READONLY)
        fgs.AddMany([(self._ctrl_log, 1, wx.EXPAND)])#, (descr, 1, wx.EXPAND)])
        self._core.Add(fgs, flag=wx.ALL|wx.EXPAND)
        self._panel.SetSizer(self._core)


    def _box_buttons(self):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        rerun_button = wx.Button(self, label='Re-run')
        rerun_button.Bind(wx.EVT_BUTTON, self._on_rerun)
        rerun_button.SetDefault()
        buttons.Add(rerun_button, 0, wx.RIGHT, border=5)
        save_button = wx.Button(self, label="Save to file")
        save_button.Bind(wx.EVT_BUTTON, self._on_save)
        buttons.Add(save_button, 0, wx.RIGHT, border=5)
        close_button = wx.Button(self, label='Close')
        close_button.Bind(wx.EVT_BUTTON, self._on_cancel)
        buttons.Add(close_button)
        self._bottom.Add(self._panel, 0, wx.EXPAND|wx.ALL, border=10)
        self._bottom.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        self._bottom.SetSizeHints(self)


    def _on_apply(self, e=None, refresh=True):
        pass


    def _on_rerun(self, e=None, refresh=True):
        log = self._ctrl_log.GetValue()
        log = log.replace('“', '"')
        log = log.replace('”', '"')
        log = log.replace('—', '--')
        #log_bck = dc(log)
        """
        i = self._gui._sess_list.index(self._gui._sess_sel)
        self._gui._panel_sess._tab.DeleteItem(i)
        del self._gui._sess_list[i]
        del self._gui._sess_item_list[i]

        # Run selected log
        self._gui._log_run(json.loads(log), skip_tab=True)
        """
        self._gui._log_rerun(log)
        self._ctrl_log.SetValue(self._log)

    def _on_cancel(self, e=None):
        self._shown = False
        self.Destroy()


    def _on_save(self, e=None, path=None):
        if path is None:
            if hasattr(self._gui, '_path'):
                path=os.path.basename(self._gui._path)
            else:
                path='.'
        name = self._gui._sess_sel.name
        with wx.FileDialog(self._gui._panel_sess, "Save log", path, name,
                           wildcard="JSON file (*.json)|*.json",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) \
                           as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            path = fileDialog.GetPath()
            dir = fileDialog.GetDirectory()
            logging.info("I'm saving log %s..." % path)
            self._gui._sess_sel.log.save(path)


    def _refresh(self):
        if self._sel != self._gui._panel_sess._sel:
            self._sel = self._gui._panel_sess._sel
        self._log = self._gui._sess_sel.log.str
        self._ctrl_log.SetValue(self._log)


class GUIDialogMiniMeta(GUIDialogMini):
    def __init__(self,
                 gui,
                 title):
        self._gui = gui
        self._gui._dlg_mini_meta = self
        self._sel = dc(self._gui._panel_sess._sel)

        self._meta = meta_parse(self._gui._sess_sel.spec._meta)
        self._meta_backup = dc(self._gui._sess_sel.spec._meta_backup)
        super(GUIDialogMiniMeta, self).__init__(gui, title)
        self.Bind(wx.EVT_CLOSE, self._on_cancel)
        self._shown = False


    def _box_ctrl(self):
        fgs = wx.FlexGridSizer(2, 1, 4, 15)
        """
        descr = wx.StaticText(
                    self._panel, -1,
                    label="Each line define a graph element as a set of\n"
                          "comma-separated values. Values are: session,\n"
                          "table, x column, y column, mask column (if any),\n"
                          "type of graph (plot, step, scatter), line style\n"
                          "or marker symbol, line width or marker size,\n"
                          "color, alpha transparency.")
        """
        self._ctrl_meta = wx.TextCtrl(self._panel, -1, value=self._meta,
                                      size=(400, 300), style = wx.TE_MULTILINE)
        #self._ctrl_z = wx.TextCtrl(self._panel, -1, value="%3.7f" % 10, size=(150, -1))
        fgs.AddMany([(self._ctrl_meta, 1, wx.EXPAND)])#, (descr, 1, wx.EXPAND)])
        self._core.Add(fgs, flag=wx.ALL|wx.EXPAND)
        self._panel.SetSizer(self._core)


    def _box_buttons(self):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        apply_button = wx.Button(self, label='Apply')
        apply_button.Bind(wx.EVT_BUTTON, self._on_apply)
        #apply_button.SetDefault()
        buttons.Add(apply_button, 0, wx.RIGHT, border=5)
        orig_button = wx.Button(self, label="Back to original")
        orig_button.Bind(wx.EVT_BUTTON, self._on_original)
        buttons.Add(orig_button)
        self._bottom.Add(self._panel, 0, wx.EXPAND|wx.ALL, border=10)
        self._bottom.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        self._bottom.SetSizeHints(self)


    def _on_apply(self, e=None, refresh=True, log=True):
        self._meta = self._ctrl_meta.GetValue()
        self._set(self._meta)
        self._gui._sess_sel.spec._meta = dc(self._meta_backup)
        #self._gui._meta_list[self._sel] = self._meta
        for m in self._meta.split('\n'):
            k = m.split(': ')[0]
            #print(k)
            v = m.split(': ')[1].split(' / ')
            self._gui._sess_sel.spec._meta[k] = v[0]
            self._gui._sess_sel.spec._meta.comments[k] = v[1]
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_meta', '_on_apply',
                                 {'e': None, 'refresh': refresh})
        if refresh: self._gui._refresh(init_cursor=True, init_tab=False)


    def _on_original(self, e=None, refresh=True, log=True):
        self._meta = meta_parse(self._meta_backup)
        self._set(self._meta)
        self._gui._sess_sel.spec._meta = dc(self._meta_backup)
        #print(meta_parse(self._gui._sess_sel.spec.meta))
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_meta', '_on_original',
                                 {'e': None, 'refresh': refresh})
        if refresh: self._gui._refresh(init_cursor=True, init_tab=False)



    def _on_cancel(self, e=None):
        self._shown = False
        self.Destroy()


    def _refresh(self):
        if self._sel != self._gui._panel_sess._sel:
            self._sel = self._gui._panel_sess._sel
            #self._elem = elem_expand(graph_elem, self._sel)
            self._meta = meta_parse(self._gui._sess_sel.spec.meta)
        self._ctrl_meta.SetValue(self._meta)
        #self._on_apply(refresh=False)


    def _set(self, value, log=True):
        if log:
            sess = self._gui._sess_sel
            sess.log.append_full('_dlg_mini_meta', '_set',
                                 {'value': value, 'log': False})
        self._ctrl_meta.SetValue(value)


class GUIDialogMiniSystems(GUIDialogMini):
    def __init__(self,
                 gui,
                 title,
                 targ=None,
                 series='CIV',
                 z=None,
                 hwin=hwin_def):
        self._gui = gui
        self._gui._dlg_mini_systems = self
        self._targ = targ
        self._series = series
        if z is None:
            try:
                z = np.mean(self._gui._sess_sel.spec.t['x'].to(au.nm))\
                        /xem_d['CIV_1548']-1
            except:
                z = 2.0
        self._z = z
        self._hwin = hwin
        super(GUIDialogMiniSystems, self).__init__(gui, title)
        self._shown = False

    def _box_ctrl(self):
        fgs = wx.FlexGridSizer(3, 2, 4, 15)
        stat_series = wx.StaticText(self._panel, -1, label="Series:")
        stat_z = wx.StaticText(self._panel, -1, label="Redshift:")
        stat_hwin = wx.StaticText(self._panel, -1, label="Half window (km/s):")
        self._ctrl_series = wx.TextCtrl(self._panel, -1, value=self._series, size=(150, -1))
        self._ctrl_z = wx.TextCtrl(self._panel, -1, value="%3.7f" % self._z, size=(150, -1))
        self._ctrl_hwin = wx.TextCtrl(self._panel, -1, value="%3.1f" % self._hwin, size=(150, -1))
        fgs.AddMany([(stat_series, 1, wx.EXPAND), (self._ctrl_series, 1, wx.EXPAND),
                     (stat_z, 1, wx.EXPAND), (self._ctrl_z, 1, wx.EXPAND),
                     (stat_hwin, 1, wx.EXPAND), (self._ctrl_hwin, 1, wx.EXPAND)])
        self._gui._sess_sel._series_sel = self._series
        self._gui._sess_sel._hwin_sel = self._hwin
        self._core.Add(fgs, flag=wx.ALL|wx.EXPAND)
        self._panel.SetSizer(self._core)

    def _box_buttons(self):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        apply_button = wx.Button(self, label='Apply')
        apply_button.Bind(wx.EVT_BUTTON, self._on_apply)
        apply_button.SetDefault()
        buttons.Add(apply_button, 0, wx.RIGHT, border=5)
        self._cursor_button = wx.Button(self, label="Show cursor")
        self._cursor_button.Bind(wx.EVT_BUTTON, self._on_show)
        buttons.Add(self._cursor_button, 0, wx.RIGHT, border=5)
        stick_button = wx.Button(self, label="Stick cursor")
        stick_button.Bind(wx.EVT_BUTTON, self._on_stick)
        buttons.Add(stick_button)
        self._bottom.Add(self._panel, 0, wx.EXPAND|wx.ALL, border=10)
        self._bottom.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        self._bottom.SetSizeHints(self)


    def _cursor_refresh(self):
        # Unfreeze cursors in case they were frozen
        self._gui._graph_main._graph._cursor_frozen = False
        if hasattr(self._gui, '_graph_det'):
            self._gui._graph_det._graph._cursor_frozen = False

        # Refresh cursor
        shown = True if self._shown else False
        self._shown = False
        self._on_cancel(e=None)
        self._shown = True
        self._on_apply(e=None)
        self._shown = shown


    def _on_apply(self, e, refresh=True):
        series = self._ctrl_series.GetValue()
        z = self._ctrl_z.GetValue()
        hwin = self._ctrl_hwin.GetValue()
        self._gui._sess_sel._series_sel = series
        self._gui._sess_sel._z_sel = float(z)
        self._gui._sess_sel._hwin_sel = float(hwin)
        if self._targ != None:
            self._targ(self._gui._sess_sel)
        if hasattr(self._gui, '_graph_det'):
            series = trans_parse(self._gui._sess_sel._series_sel)
            self._gui._graph_det._graph._fig.clear()
            self._gui._graph_det._update(series, float(z), float(hwin))
        if refresh: self._gui._refresh(init_cursor=True, init_tab=False)


    def _on_show(self, e):
        sel = self._gui._graph_main._sel
        if not self._shown:
            sel.append(self._gui._cursor.key)
            self._on_apply(e)
            self._cursor_button.SetLabel("Hide cursor")
        else:
            sel.remove(self._gui._cursor.key)
            self._on_cancel(e)
            self._cursor_button.SetLabel("Show cursor")
        self._shown = not self._shown


    def _on_stick(self, e):
        self._gui._graph_main._on_cursor_stick(e)


    def _on_cancel(self, e):
        if hasattr(self._gui, '_cursor'):
            self._gui._cursor.Check(False)
            if hasattr(self._gui, '_graph_det') \
                and hasattr(self._gui._graph_det._graph, '_cursor'):
                del self._gui._graph_det._graph._cursor
        if hasattr(self._gui._sess_sel, '_series_sel'):
            del self._gui._sess_sel._series_sel
        if hasattr(self._gui._sess_sel, '_hwin_sel'):
            del self._gui._sess_sel._hwin_sel
        self._gui._refresh(init_cursor=True, init_tab=False)


    def _refresh(self, series='CIV', z=2.0):
        self._ctrl_series.SetValue(series)
        self._ctrl_z.SetValue("%3.7f" % z)
        if hasattr(self._gui._graph_det._graph, '_cursor'):
            self._on_apply(None)
