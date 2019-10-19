from collections import OrderedDict
from copy import deepcopy as dc
import inspect
import numpy as np
import wx

class GUIDialog(wx.Dialog):

    def __init__(self,
                 gui,
                 title):
        self._gui = gui
        self._gui._dlg = self

        super(GUIDialog, self).__init__(parent=None, title=title)

class GUIDialogMethod(wx.Dialog):

    def __init__(self,
                 gui,
                 title,
                 attr):

        self._gui = gui
        self._gui._dlg_method = self
        super(GUIDialogMethod, self).__init__(parent=None, title=title)
        self._attr = np.array(attr, ndmin=1)
        self._methods = []
        self._params = []
        self._brief = []
        self._details = []
        self._doc = []

        for a in self._attr:
            obj = self._gui._sess_sel
            method = getattr(obj, a)
            self._methods.append(method)
            self._get_params(method)
            self._get_doc(method)
        self._init()

    def _init(self, attr=None, methods=None, params=None, brief=None,
              details=None, doc=None):
        if attr == None: attr = self._attr
        if methods == None: methods = self._methods
        if params == None: params = self._params
        if brief == None: brief = self._brief
        if details == None: details = self._details
        if doc == None: doc = self._doc

        print(attr, methods, params, brief, details, doc)
        self._ctrl = []
        panel = wx.Panel(self)
        bottom = wx.BoxSizer(wx.VERTICAL)
        core = wx.BoxSizer(wx.VERTICAL)
        if len(methods) == 1:
            self._box_descr(panel, core, details)
            self._box_params(panel, core, params, doc)
        else:
            self._box_methods(panel, core, attr, methods, params, brief,
                              details, doc)
        self._box_buttons(panel, bottom)
        self.SetSizer(bottom)
        self.Centre()
        self.Show()

    def _box_buttons(self, panel, box):
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        cancel_button = wx.Button(self, label='Cancel')
        cancel_button.Bind(wx.EVT_BUTTON, self._on_cancel)
        run_button = wx.Button(self, label='Run')
        run_button.Bind(wx.EVT_BUTTON, self._on_run)
        run_button.SetDefault()
        buttons.Add(cancel_button, 0, wx.RIGHT, border=5)
        buttons.Add(run_button, 0)

        box.Add(panel, 0, wx.EXPAND|wx.ALL, border=10)
        box.Add(buttons, 0, wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT|wx.BOTTOM,
                     border=10)
        box.SetSizeHints(self)



    def _box_descr(self, panel, box, details):

        # Description
        sb = wx.StaticBox(panel, label="Description")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        st = wx.StaticText(sb, 1, label='\n'.join(details))
        st.Wrap(400)
        sbs.Add(st, flag=wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND, border=8)
        box.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        panel.SetSizer(box)

    def _box_methods(self, panel, box, attr, methods, params, brief, details,
                     doc):
        sb = wx.StaticBox(panel, label="Methods")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        bbs = wx.BoxSizer(wx.VERTICAL)
        for a_l, m_l, p_l, b_l, de_l, d_l \
            in zip(attr, methods, params, brief, details, doc):
            button = wx.Button(panel, label=b_l)
            # Both method and attribute should be argument of lambda, otherwise
            # only the last method and attribute are passed
            button.Bind(
                wx.EVT_BUTTON,
                lambda e, a=[a_l], m=[m_l], p=[p_l], b=[b_l], de=[de_l], d=[d_l]: \
                self._on_method(e, a, m, p, b, de, d))
            bbs.Add(button, 0, wx.BOTTOM|wx.EXPAND, border=5)
        sbs.Add(bbs, flag=wx.LEFT|wx.TOP, border=20)
        box.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        panel.SetSizer(box)

    def _box_params(self, panel, box, params, doc):
        # Parameters
        sb = wx.StaticBox(panel, label="Parameters")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        len_params = np.sum([len(i) for i in params])
        fgs = wx.FlexGridSizer(len_params, 2, 4, 15)
        fgs_add = []
        for p_l, d_l in zip(params, doc):
            ctrl_l = []
            for p, d in zip(p_l, d_l):
                stat = wx.StaticText(panel, -1, label=d+':')
                ctrl = wx.TextCtrl(panel, -1, value=str(p_l[p]))
                fgs_add.append((stat, 1, wx.EXPAND))
                fgs_add.append((ctrl, 1, wx.EXPAND))
                ctrl_l.append(ctrl)
            self._ctrl.append(ctrl_l)
        if np.size(self._ctrl) > 0:
            fgs.AddMany(fgs_add)
            sbs.Add(fgs, flag=wx.ALL|wx.EXPAND, border=8)
            box.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        panel.SetSizer(box)


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

    def _on_method(self, e, attr, method, params, brief, details, doc):
        self._init(attr, method, params, brief, details, doc)

    def _on_run(self, e):
        for a, p_l, c_l in zip(self._attr, self._params, self._ctrl):
            obj = self._gui._sess_sel
            m = getattr(obj, a)

            for p, c in zip(p_l, c_l):
                pmod = c.GetValue()
                p_l[p] = pmod
            out = m(**p_l)
            if out is not None:
                if out is 0:
                    self._gui._panel_sess._refresh()
                    self._gui._graph_main._refresh(self._gui._sess_items)
                else:
                    self._gui._panel_sess._on_add(out, open=False)
                self.Close()
