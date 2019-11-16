from collections import OrderedDict
from copy import deepcopy as dc
import inspect
import numpy as np
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
        for a in self._attr:
            if self._obj == None:
                self._obj = self._gui._sess_sel
            method = getattr(self._obj, a)
            self._methods.append(method)
            self._get_params(method)
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

    def _get_params(self, method):
        keys = inspect.getargspec(method)[0][1:]
        defs = inspect.getargspec(method)[-1]
        if defs == None:
            defs = []
        defs = [str(d) for d in defs]
        values = np.append(['']*(len(keys)-len(defs)), defs)
        self._params.append(OrderedDict(zip(keys, values)))

    def _on_cancel(self, e):
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

            out = m(**p_l)
            if out is not None:
                if out is 0:
                    self._gui._panel_sess._refresh()
                    self._gui._panel_sess._menu._refresh()
                    self._gui._graph_main._refresh(self._gui._sess_items)
                    if hasattr(self._gui, '_graph_det'):
                        xlim = self._gui._graph_det._graph._ax.get_xlim()
                        ylim = self._gui._graph_det._graph._ax.get_ylim()
                        self._gui._graph_det._refresh(self._gui._sess_items, xlim=xlim,
                                                      ylim=ylim)
                else:
                    self._gui._panel_sess._on_add(out, open=False)
                self.Close()

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
                 params_parent=None):

        super(GUIDialogMethod, self).__init__(gui, title, attr, obj)
        self._cancel_run = cancel_run
        self._params_parent = params_parent
        self._box_descr()
        self._box_params()
        self._box_buttons(self._cancel_run)
        self.SetSizer(self._bottom)
        self.Centre()
        self.SetPosition((self.GetPosition()[0], wx.DisplaySize()[1]*0.25))
        self.Show()

    def _box_descr(self):

        # Description
        sb = wx.StaticBox(self._panel, label="Description")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        st = wx.StaticText(sb, 1, label='\n'.join(self._details))
        st.Wrap(400)
        sbs.Add(st, flag=wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND, border=8)
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
        self._panel.SetSizer(self._core)


class GUIDialogMethods(GUIDialog):

    def __init__(self,
                 gui,
                 title,
                 attr,
                 obj=None):

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
                        params_parent=self._params)
