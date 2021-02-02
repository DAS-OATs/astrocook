from .functions import elem_expand
from .graph import Graph
from .gui_dialog import * #GUIDialogMini
from .syst_list import SystList
from .vars import *
from collections import OrderedDict
import logging
from matplotlib import pyplot as plt
from matplotlib.backend_bases import Event
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
from matplotlib.figure import Figure
import numpy as np
from scipy.stats import norm
import wx

class GUIGraphMain(wx.Frame):
    """ Class for the GUI spectrum graph frame """

    def __init__(self,
                 gui,
                 #sess,
                 title="Spectrum",
                 size_x=wx.DisplaySize()[0]*0.96,
                 size_y=wx.DisplaySize()[1]*0.5,
                 main=True,
                 **kwargs):
        """ Constructor """

        self._gui = gui
        self._title = title
        self._size_x = size_x
        self._size_y = size_y
        self._sel = graph_sel
        self._cols_sel = graph_cols_sel
        sel = len(self._gui._sess_list)
        if sel>0: sel -= 1
        self._elem = elem_expand(graph_elem, sel)

        self._logx = False
        self._logy = False
        self._norm = False
        self._legend = False
        self._closed = False
        if main:
            self._gui._graph_main = self
        self._init(**kwargs)
        self.SetPosition((wx.DisplaySize()[0]*0.02, wx.DisplaySize()[1]*0.40))


    def _init(self, **kwargs):
        super(GUIGraphMain, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))

        self._panel = wx.Panel(self)
        self._graph = Graph(self._panel, self._gui, self._sel, **kwargs)
        self._textbar = wx.StaticText(self._panel, 1, style=wx.ALIGN_RIGHT)
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        box_toolbar.Add(self._graph._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._textbar, 1, wx.RIGHT)
        self._box = wx.BoxSizer(wx.VERTICAL)
        self._box.Add(self._graph._canvas, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP)
        self._panel.SetSizer(self._box)
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)
        #self._graph._toolbar.wx_ids['Home'] = self._home

        #self._gui._statusbar = self.CreateStatusBar()
        move_id = self._graph._canvas.mpl_connect('motion_notify_event',
                                                  self._graph._on_move)
        click_id = self._graph._canvas.mpl_connect('button_release_event',
                                                   self._graph._on_click)


    def _refresh(self, sess, **kwargs):
        if self._closed:
            self._init()

        self._graph._refresh(sess, self._logx, self._logy, self._norm,
                             self._legend, **kwargs)
        self.Show()

    #def _on_line_new(self, event):
    #    print(self._click_xy)


    def _on_cursor_stick(self, event=None, cursor_z=None):
        sess = self._gui._sess_sel
        #print(self._graph._cursor._z * sess.spec._rfz)
        z = "%2.6f" % (self._graph._cursor._z*(1+sess.spec._rfz))

        if not hasattr(sess, '_cursors'):
            sess._cursors = {z: self._graph._cursor}
        else:
            sess._cursors[z] = self._graph._cursor
        sess._graph_elem += \
            '\n%i,cursor,%s,None,None,axvline,:,1.0,C%s,1.0' \
            % (self._gui._panel_sess._sel, z, (len(sess._cursors)-1)%10)
        self._elem = sess._graph_elem

        if hasattr(self._gui, '_dlg_mini_graph'):
            self._gui._dlg_mini_graph._refresh()
        self._refresh(sess)


    def _on_node_add(self, event):
        sess = self._gui._sess_sel
        x, y = sess._clicks[-1][0], sess._clicks[-1][1]
        sess.log.append_full('cb', 'node_add', {'x': x, 'y': y})
        sess.log.append_full('cb', 'nodes_interp', {})
        sess.spec._node_add(sess.nodes, x, y)
        sess.spec._nodes_interp(sess.lines, sess.nodes)
        self._gui._refresh()

    def _on_node_remove(self, event):
        sess = self._gui._sess_sel
        x, y = sess._clicks[-1][0], sess._clicks[-1][1]
        sess.log.append_full('cb', 'node_remove', {'x': x})
        sess.log.append_full('cb', 'nodes_interp', {})
        sess.spec._node_remove(sess.nodes, x)
        sess.spec._nodes_interp(sess.lines, sess.nodes)
        self._gui._refresh()

    def _on_region_extract(self, event):
        sess = self._gui._sess_sel
        x = [sess._clicks[0][0], sess._clicks[1][0]]
        xmin = np.min(x)
        xmax = np.max(x)
        sel_old = self._gui._sess_list.index(sess)
        reg = sess.cb.region_extract(xmin, xmax)
        #self._gui._refresh()
        self._gui._panel_sess._on_add(reg, open=False)

        sess = self._gui._sess_sel
        sess_list = [self._gui._sess_list[sel_old]]
        sess.log.merge_full('cb', 'region_extract',
                             {'xmin': xmin, 'xmax': xmax}, sess_list, sess)

    def _on_spec_zap(self, event):
        sess = self._gui._sess_sel
        x = [sess._clicks[0][0], sess._clicks[1][0]]
        xmin = np.min(x)
        xmax = np.max(x)
        sess.log.append_full('cb', 'feature_zap', {'xmin': xmin, 'xmax': xmax})
        sess.spec._zap(xmin, xmax)
        self._gui._refresh()

    def _on_stats_show(self, event):
        sess = self._gui._sess_sel
        try:
            x = [sess._clicks[-2][0], sess._clicks[-1][0]]
        except:
            x = (0, np.inf)
            sess._shade = False
        xmin = np.min(x)
        xmax = np.max(x)
        sess.spec._stats_print(xmin, xmax)
        sess._clicks = []
        sess._stats = True
        #self._graph._stats = True
        self._gui._refresh()

    def _on_stats_hide(self, event):
        sess = self._gui._sess_sel
        del sess.spec._stats_text_red
        sess._clicks = []
        sess._stats = False
        sess._shade = False
        self._gui._refresh()

    def _on_syst_new(self, event):
        sess = self._gui._sess_sel

        params = [{'series': sess._series_sel, 'z': self._graph._cursor._z,
                   'refit_n': 0, 'sel_fit': True}]
        dlg = GUIDialogMethod(self._gui, 'New system', 'syst_new',
                              params_last = params)
        self._gui._refresh(init_cursor=True)


    def _on_close(self, event=None):
        self._closed = True
        self.Destroy()
        del self._gui._graph_det


class GUIGraphDetail(GUIGraphMain):

    def __init__(self,
                 gui,
                 #sess,
                 title="Spectrum detail",
                 size_x=wx.DisplaySize()[0]*0.4,
                 size_y=wx.DisplaySize()[1]*0.4,
                 **kwargs):
        super(GUIGraphDetail, self).__init__(gui, title, size_x, size_y,
                                             main=False, **kwargs)
        self._norm = True
        self._gui._graph_det = self
        #self._gui._sess_sel._graph_det = self
        self._graph._legend = False
        self._graph._cursor_lines = []
        sel = len(self._gui._sess_list)
        if sel>0: sel -= 1
        self._elem = elem_expand(graph_elem, sel)
        #self.SetPosition((wx.DisplaySize()[0]*0.58, wx.DisplaySize()[0]*0.02))

    def _define_lim(self, x, t=None, xspan=30, ymargin=0.1):#, norm=False):
        if t == None:
            t = self._gui._sess_sel.spec.t
            t = t[np.logical_and(~np.isnan(t['x']), ~np.isnan(t['y']))]
        w = np.argmin(np.abs(t['x']-x))
        xmin = t['x'][max(w-xspan, 0)]
        xmax = t['x'][min(w+xspan, len(t)-1)]
        ysel = t['y'][np.where(np.logical_and(t['x']>xmin, t['x']<xmax))]
        yspan = ymargin*(np.max(ysel)-np.min(ysel))
        xlim = (xmin, xmax)
        #if norm:
        ylim = (-0.2, 1.2)
        #else:
        #    ylim = (np.min(ysel)-yspan, np.max(ysel)+yspan)

        return xlim, ylim

    def _update(self, series, z, hwin, nmax=15):
        graph = self._graph
        if len(series) > nmax:
            logging.warning("I'm discarding the last %i transitions, because I "
                            "can display only 15 panels." % (len(series)-nmax))
            series = series[:nmax]
        n = len(series)
        rows = min(5, n)
        cols = (n-1)//5+1
        idxs = np.ravel(np.reshape(range(rows*cols), (rows, cols)).T)[:n]
        size_x = min(wx.DisplaySize()[0]*0.5, wx.DisplaySize()[0]*0.4*cols)
        size_y = min(wx.DisplaySize()[1]*0.9, wx.DisplaySize()[1]*0.3*rows)
        self.SetSize(wx.Size(size_x, size_y))
        self.SetPosition((min(wx.DisplaySize()[0]*0.98-size_x,
                              wx.DisplaySize()[0]*0.58),
                         wx.DisplaySize()[1]*0.02))

        # Redshift and wavelengths need to be initialized before the cursor
        # is created in the graph
        graph._axes = OrderedDict()
        graph._zems = OrderedDict()
        graph._xs = OrderedDict()
        graph._series = OrderedDict()
        #graph._cursor_lines = []
        for i, s in zip(idxs, series):
            #key = s[-4:]
            key = s.split('_')[-1]
            key = s
            x = (1+z)*xem_d[s]
            zem = (1+z)*xem_d[s]/xem_d['Ly_a']-1
            #print(x)
            #print(zem)
            #print('out', xem_d[s], graph._zem, graph._x)
            graph._zems[key] = zem
            graph._xs[key] = x
            graph._series[key] = s
            graph._z = z
            if i == 0:
                graph._x = x
                graph._zem = zem

        #graph._axes = []
        #for i, (x, zem) in enumerate(zip(graph._xs, graph._zems)):
            xunit = self._gui._sess_sel.spec.x.unit
            self._gui._sess_sel.cb.x_convert(zem=zem)
            self._gui._sess_sel._xdet = x
            self._gui._sess_sel._ydet = 0.0
            _, ylim = self._define_lim(0)#, norm=True)

            if i == 0:
                graph._ax = graph._fig.add_subplot(rows, cols, i+1)
                #title = series
                #if len(series) > 1:
                #    graph._ax.tick_params(labelbottom=False)
            else:
                title = None
                graph._ax = graph._fig.add_subplot(rows, cols, i+1,
                                                   sharex=graph._ax,
                                                   sharey=graph._ax)
            if i < len(series)-cols:
                graph._ax.tick_params(labelbottom=False)
            if i%cols !=0:
                graph._ax.tick_params(labelleft=False)

            graph._ax.tick_params(top=True, right=True, direction='in')
            graph._fig.subplots_adjust(hspace=0.,wspace=0.)
            graph._axes[key] = graph._ax
            self._refresh(
                self._gui._sess_items, text=key,
                xlim=(-hwin,hwin), ylim=ylim)

            self._gui._sess_sel.cb.x_convert(zem=zem, xunit=xunit)


class GUIGraphHistogram(GUIGraphMain):

    def __init__(self,
                 gui,
                 title="Column histogram",
                 size_x=wx.DisplaySize()[0]*0.6,
                 size_y=wx.DisplaySize()[1]*0.6,
                 **kwargs):
        #self._col_values = col_values
        super(GUIGraphHistogram, self).__init__(gui, title, size_x, size_y,
                                                main=False, **kwargs)
        self._gui._graph_hist = self
        self.SetPosition((wx.DisplaySize()[0]*0.48, wx.DisplaySize()[1]*0.02))

    def _init(self, **kwargs):
        super(GUIGraphMain, self).__init__(parent=None, title=self._title,
                                           size=(self._size_x, self._size_y))


        self._panel = wx.Panel(self)
        #self._graph = Graph(self._panel, self._gui, self._sel, **kwargs)
        self._fig = Figure()
        self._canvas = FigureCanvasWxAgg(self._panel, -1, self._fig)
        self._toolbar = NavigationToolbar2Wx(self._canvas)
        self._toolbar.Realize()
        self._textbar = wx.StaticText(self._panel, 1, style=wx.ALIGN_RIGHT)
        box_toolbar = wx.BoxSizer(wx.HORIZONTAL)
        #box_toolbar.Add(self._graph._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._toolbar, 1, wx.RIGHT)
        box_toolbar.Add(self._textbar, 1, wx.RIGHT)
        self._box = wx.BoxSizer(wx.VERTICAL)
        #self._box.Add(self._graph._canvas, 1, wx.EXPAND)
        self._box.Add(self._canvas, 1, wx.EXPAND)
        self._box.Add(box_toolbar, 0, wx.TOP)
        self._panel.SetSizer(self._box)
        self.Centre()
        self.Bind(wx.EVT_CLOSE, self._on_close)

    def _refresh(self, sess, **kwargs):
        if self._closed:
            self._init()
        """
        self._graph._ax.clear()
        self._graph._init_ax(111)
        self._graph._ax.hist(col_values)
        self._graph._canvas.draw()
        """
        if hasattr(self, '_ax'):
            self._ax.clear()
        self._ax = self._fig.add_subplot(111)
        self._ax.tick_params(top=True, right=True, direction='in')#self._init_ax(111)
        label = self._gui._col_tab.GetColLabelValue(self._gui._col_sel).split('\n')[0]
        self._ax.set_xlabel(label.replace('\n',' ')+' km/s')
        self._ax.set_ylabel('Frequency')
        values = self._gui._col_values

        """
        # Temporary: only for deltav
        values = []
        rej = 0
        for v in self._gui._col_values:
            if np.abs(v) < 19.98 and v != 0 and v not in values:
                values.append(v)
            else:
                rej += 1
        """
        """
        scale = int(round(np.median(np.log10(values))))-1
        min = np.floor(np.min(values))
        max = np.ceil(np.max(values))
        if np.log10(max)-np.log10(min) > 2:
            values = np.log10(values)
            scale = int(round(np.median(np.log10(values))))-1
            min = np.floor(np.min(values))
            max = np.ceil(np.max(values))
        step = 0.5*10**scale
        bins = np.arange(min-0.5*step, max+1.5*step, step)
        """
        scale = np.ceil(np.log(np.abs((np.max(values)-np.min(values))/np.median(values))))
        bins = int(scale)*10
        step=2/3
        #min = -5
        min = np.floor(np.min(values))
        #max = 5
        max = np.ceil(np.max(values))
        bins = np.arange(min-0.0*step, max+1.0*step, step)
        n, bins, patches = self._ax.hist(values, align='mid', bins=bins)
        #print(n, bins)
        #mu = np.average(bins[:-1]+0.25, weights=n)
        #sigma = np.sqrt(np.average((bins[:-1]+0.25-mu)**2, weights=n))
        mu = np.mean(values)
        sigma = np.std(values)
        sigmal = np.percentile(values, 15.87)
        sigmar = np.percentile(values, 84.13)
        clr = 0.5*(sigmal+sigmar)
        dv = 0.83
        x = np.linspace(bins[0], bins[-1], len(bins)*10)
        #g = np.exp(-0.5 * ((x-clr) / dv)**2)
        g = np.zeros(x.shape)
        if 'sigmav_AB' in self._gui._sess_sel.systs.t.colnames:
            col = 'sigmav'
        elif 'sigmav' in self._gui._sess_sel.systs.t.colnames:
            col = 'sigmav'
        else:
            col = None
        if col != None:
            for sigmav in self._gui._sess_sel.systs.t[col]:
                if not np.isnan(sigmav):
                    g = g+np.exp(-0.5 * ((x-clr) / sigmav)**2)
        g = g*len(g)/np.sum(g) * np.sum(n)/len(bins)
        #print(np.sum(n)*(bins[1]-bins[0]), np.sum(g)/len(g))
        self._ax.axvline(mu, linestyle='--', c='C1', label=r'$\mu$ = %3.2f, $\sigma$ = %3.2f'%(mu,sigma))
        self._ax.axvline(mu+sigma, linestyle='--', c='C1')
        self._ax.axvline(mu-sigma, linestyle='--', c='C1')
        self._ax.axvline(sigmal, linestyle=':', c='C2', label='15.87th cent = %3.2f,\n84.13th cent = %3.2f'%(sigmal,sigmar))
        self._ax.axvline(sigmar, linestyle=':', c='C2')
        #self._ax.plot(x, g, c='C3', label=r'$\langle\Delta v\rangle_{\mathrm{SNR}=10} = %3.2f$'%dv)
        #self._ax.plot(x, g, c='C3', label=r'$\sqrt{\langle\Delta v\rangle_{\mathrm{SNR}=10}^2+\langle\Delta v\rangle_{\mathrm{SNR}=15}^2} = %3.2f$'%dv)
        self._ax.plot(x, g, c='C3', label="Cumulative position uncertainty")
        """
        self._ax.text(0.95, 0.9, r'$\mu$ = %3.2f' % mu, c='C1',
                      transform=self._ax.transAxes, horizontalalignment='right')
        self._ax.text(0.95, 0.8, r'$\sigma$ = %3.4f' % sigma, c='C1',
                      transform=self._ax.transAxes, horizontalalignment='right')
        self._ax.text(0.95, 0.7, '15.87th cent = %3.4f' % sigmal, c='C2',
                      transform=self._ax.transAxes, horizontalalignment='right')
        self._ax.text(0.95, 0.6, '84.13th cent = %3.4f' % sigmar, c='C2',
                      transform=self._ax.transAxes, horizontalalignment='right')
        self._ax.text(0.95, 0.5, r'$\sqrt{\langle\Delta v\rangle_{\mathrm{SNR}=10}^2+\langle\Delta v\rangle_{\mathrm{SNR}=15}^2} = %3.2f$' % dv, c='C3',
                      transform=self._ax.transAxes, horizontalalignment='right')
        """
        self._ax.legend()
        self._canvas.draw()
        self.Show()
