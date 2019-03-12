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
                 source,
                 targ,
                 attr,
                 cl=None):

        self._gui = gui
        self._gui._dlg_method = self
        super(GUIDialogMethod, self).__init__(parent=None, title=title)
        self._source = np.array(source, ndmin=1)
        self._targ = np.array(targ, ndmin=1)
        self._attr = np.array(attr, ndmin=1)
        if cl != None:
            self._cl = np.array(cl, ndmin=1)
        else:
            self._cl = np.array([None]*len(self._source))
        self._methods = []
        self._params = []
        self._brief = []
        self._doc = []

#        if cl is not None:
#            setattr(self._gui._sess_sel, source, cl(sess=self._gui._sess_sel))


        for s, a, c in zip(self._source, self._attr, self._cl):
            if c != None and getattr(self._gui._sess_sel, s) == None:
                print(s, c)
                setattr(self._gui._sess_sel, s, c(sess=self._gui._sess_sel))

            if s == None:
                obj = self._gui._sess_sel
            else:
                obj = getattr(self._gui._sess_sel, s)
            method = getattr(obj, a)
            self._methods.append(method)
            #super(GUIDialogMethod, self).__init__(parent=None, title=title)
            self._get_params(method)
            self._get_doc(method)


        panel = wx.Panel(self)
        box = wx.BoxSizer(wx.VERTICAL)
        core = wx.BoxSizer(wx.VERTICAL)

        # Description
        sb = wx.StaticBox(panel, label="Description")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        st = wx.StaticText(sb, 1, label='\n'.join(self._brief))
        st.Wrap(400)
        sbs.Add(st, flag=wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND, border=8)
        core.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        panel.SetSizer(core)


        # Parameters
        sb = wx.StaticBox(panel, label="Parameters")
        sbs = wx.StaticBoxSizer(sb, wx.VERTICAL)
        len_params = np.sum([len(i) for i in self._params])
        fgs = wx.FlexGridSizer(len_params, 2, 4, 15)
        fgs_add = []
        self._ctrl = []
        for p_l, d_l in zip(self._params, self._doc):
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
            core.Add(sbs, flag=wx.ALL|wx.EXPAND, border=5)
        panel.SetSizer(core)

        # Buttons
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

        self.SetSizer(box)
        self.Centre()
        self.Show()

    def _get_doc(self, method):
        full = inspect.getdoc(method)
        split = full.split('@')
        self._brief.append([s[6:-1] for s in split if s[0:5]=='brief'][0].replace('\n', ' '))
        self._doc.append([s[6:-1].split(' ', 1)[1] \
                          for s in split if s[0:5]=='param'])

    def _get_params(self, method):
        keys = inspect.getargspec(method)[0][1:]
        defs = inspect.getargspec(method)[-1]
        if defs == None:
            defs = []
        values = np.append(['']*(len(keys)-len(defs)), defs)
        self._params.append(OrderedDict(zip(keys, values)))

    def _on_cancel(self, e):
        self.Close()

    def _on_run(self, e):
        for s, t, a, p_l, c_l in zip(self._source, self._targ, self._attr,
                                     self._params, self._ctrl):
            if s == None:
                obj = self._gui._sess_sel
            else:
                obj = getattr(self._gui._sess_sel, s)
            m = getattr(obj, a)

            for p, c in zip(p_l, c_l):
                pmod = c.GetValue()
                p_l[p] = pmod
            out = m(**p_l)
            if out is not None:
                if out is 0:
                    self._gui._panel_sess._refresh()
                    self._gui._graph_spec._refresh(self._gui._sess_items)
                else:
                    if t == None:
                        new_sess = out
                    else:
                        new_sess = dc(self._gui._sess_sel)
                        setattr(new_sess, t, out)
                    self._gui._panel_sess._on_add(new_sess, open=False)
                self.Close()
