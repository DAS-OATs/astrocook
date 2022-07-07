from .message import *
from .functions import trans_parse, x_convert
from .vars import *
from astropy import units as au
from astropy import constants as aconst
from copy import deepcopy as dc
import logging
import matplotlib
matplotlib.use('WxAgg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg, _convert_agg_to_wx_bitmap
from matplotlib.figure import Figure
import matplotlib.ticker as mticker
import matplotlib.transforms as transforms
from matplotlib.widgets import Cursor
import numpy as np
import wx

# Force a given format in axis - currently not uses
# From https://stackoverflow.com/questions/49351275/matplotlib-use-fixed-number-of-decimals-with-scientific-notation-in-tick-labels
class ScalarFormatterForceFormat(mticker.ScalarFormatter):
    def _set_format(self):  # Override function that finds format to use.
        self.format = "%1.2f"  # Give format here

class Graph(object):

    def __init__(self, panel, gui, sel, init_canvas=True, init_ax=True):
        self._panel = panel
        self._gui = gui
        self._sel = sel
        self._fig = Figure()
        self._cursor_lines = []
        self._cursor_frozen = False
        self._zoom = False
        self._click_1 = False

        if init_canvas:
            self._init_canvas()
        if init_ax:
            self._init_ax(111)

    def _init_ax(self, *args):
        self._ax = self._fig.add_subplot(*args)
        self._ax.tick_params(top=True, right=True, direction='in')
        #self._cursor = Cursor(self._ax, useblit=True, color='red',
        #                      linewidth=0.5)
        #cid =  plt.connect('motion_notify_event', self._on_move)
        self._axt = None

    def _init_canvas(self):
        #self._c = 0
        #self._fig.tight_layout()#rect=[-0.03, 0.02, 1.03, 1])
        self._canvas = FigureCanvasWxAgg(self._panel, -1, self._fig)
        self._toolbar = NavigationToolbar2Wx(self._canvas)
        #self._cursor = Cursor(self._ax, useblit=True, color='red',
        #                      linewidth=0.5)
        self._toolbar.Realize()
        if self._panel is self._gui._graph_main._panel:
            focus = self._gui._graph_main
        if hasattr(self._gui, '_graph_det'):
            if self._panel is self._gui._graph_det._panel:
                focus = self._gui._graph_det
        #cid =  plt.connect('motion_notify_event', self._on_move)

    def _check_units(self, sess, axis='x'):
        unit = axis+'unit'
        _unit = '_'+unit
        if getattr(getattr(sess, 'spec'), _unit) != getattr(self, _unit) \
            and getattr(self, _unit) != au.dimensionless_unscaled:
            logging.info("I'm converting the %s unit of %s to plot it over the "
                         "data already present." % (axis, sess.name))
            getattr(sess.cb, axis+'_convert')(**{unit: getattr(self, _unit)})
            self._gui._panel_sess._refresh()

    def _on_click(self, event):
        if not event.inaxes: return
        x = float(event.xdata)
        y = float(event.ydata)
        sess = self._gui._sess_sel
        x = x/(1+sess.spec._rfz)
        from .gui_table import GUITablePopup
        if self._panel is self._gui._graph_main._panel:
            focus = self._gui._graph_main
        if hasattr(self._gui, '_graph_det'):
            if self._panel is self._gui._graph_det._panel:
                focus = self._gui._graph_det
        focus._click_xy = (x,y)
        title = []
        attr = []

        if event.button == 1:
            sess._clicks = [(x,y)]
            self._click_1 = True
        if event.button == 3:
            if self._click_1:
                sess._clicks.append((x,y))
            else:
                sess._clicks = [(x,y)]
            sess._click_1 = False

        if event.button == 3:
            title.append('Zap bin')
            attr.append('bin_zap')
            if focus == self._gui._graph_main:
                title.append('Show stats')
                attr.append('stats_show')
            if sess._stats and focus == self._gui._graph_main:
                title.append('Hide stats')
                attr.append('stats_hide')
            if len(sess._clicks) > 1 and focus == self._gui._graph_main:
                title.append('Extract region')
                attr.append('region_extract')
                title.append('Zap feature')
                attr.append('spec_zap')
                self._reg_shade()
            if hasattr(self._gui._sess_sel, 'nodes') \
                and hasattr(self._gui._sess_sel.nodes, 'x') \
                and focus == self._gui._graph_main:
                nodes = self._gui._sess_sel.nodes
                dist_x = np.abs(nodes.x.to(nodes._xunit).value-x).min()
                dist_mean = np.mean(nodes.x[1:]-nodes.x[:-1]).to(nodes._xunit).value
                title.append('Add node')
                attr.append('node_add')
                if dist_x < 0.1*dist_mean:
                    title.append('Remove node')
                    attr.append('node_remove')
            if 'cursor_z_series' in self._sel:
                title.append('Stick cursor')
                attr.append('cursor_stick')
            if 'cont' in sess.spec._t.colnames \
                and 'cursor_z_series' in self._sel:
                title.append('New system')
                attr.append('syst_new')

        focus.PopupMenu(
            GUITablePopup(self._gui, focus, event,
                                  #['Add lines', 'Add system'],
                                  #['add_line', 'add_syst']))
                                  #'Add system', 'add_syst'))
                                  title, attr))

        """
        elif 'cursor_z_series' in self._sel:
            focus.PopupMenu(GUITablePopup(self._gui, focus,
                                          event, 'Add lines', 'add_line'))
        else:
            focus.PopupMenu(GUITablePopup(self._gui, focus,
                                          event, 'Add line', 'add_line'))
        """

    def _on_zoom(self, event):
        self._zoom = True


    def _on_move(self, event):
        if not event.inaxes: return
        x = float(event.xdata)
        y = float(event.ydata)

        # Identify which axis was clicked on; compute shift in x based for
        # detail graph with respect to the last axis
        if hasattr(self, '_axes'):
            klast = tuple(self._axes)[-1]
            for k in self._axes:
                if self._axes[k] == event.inaxes:
                    ax = self._axes[k]
                    dx = aconst.c*(xem_d[k]/xem_d[klast]-1)
        else:
            ax = self._ax
            dx = 0*au.nm


        sess = self._gui._sess_sel
        x = x/(1+sess.spec._rfz)
        if self._panel is self._gui._graph_main._panel:
            focus = self._gui._graph_main
        if hasattr(self._gui, '_graph_det'):
            if self._panel is self._gui._graph_det._panel:
                focus = self._gui._graph_det
        #print(self._cursor_lines)

        """
            for curve in self._ax.get_lines():
                print(curve)
                if curve.contains(event)[0]:
                    try:
                        self._ciao.remove()
                    except:
                        pass
                    self._ciao = self._ax.text(x,y, "ciao")
        """

        # Make system id appear when you hover close enough to a system axvline
        try:
            if self._systs_id:
                xdiff = np.abs((self._systs_x-dx.to(self._systs_x.unit)).value-x)
                argmin = np.argmin(xdiff)
                try:
                    self._tag.remove()
                except:
                    pass
                if self._systs_x.si.unit == au.m: thres = 0.5
                if self._systs_x.si.unit == au.m/au.s: thres = 5
                if xdiff[argmin] < thres:
                    self._tag = ax.text(x,y, "%s\nx = %1.7f %s\nz = %1.7f" \
                                        % (self._systs_series[argmin],
                                           self._systs_l[argmin].value,
                                           self._systs_l[argmin].unit,
                                           self._systs_z[argmin]),
                                        color=self._systs_color)
        except:
            pass

        if 'cursor_z_series' in self._sel and self._cursor_frozen == False:
            if hasattr(self, '_xs'):
                for l, key in zip(self._cursor_lines, self._xs):
                    if self._text != None:
                        z = self._z+(1+self._z)*x/aconst.c.to(au.km/au.s).value
                        self._cursor._x = (self._cursor._xem*(1+z)*au.nm).to(sess.spec._xunit)
                        xem = self._xs[key]
                        self._cursor._x = (np.log(self._cursor._x/xem))*aconst.c.to(au.km/au.s)
                    else:
                        z = x/self._cursor._xmean.to(sess.spec._xunit).value-1
                        self._cursor._x = (self._cursor._xem*(1+z)*au.nm).to(sess.spec._xunit)
                    for c, xi in zip(l, self._cursor._x):
                        c.set_xdata(xi)
                        c.set_alpha(0.5)
                    #self._canvas.draw()
                    z_obs = z*(1+sess.spec._rfz)
                    focus._textbar.SetLabel("x=%2.4f, y=%2.4e; z[%s]=%2.5f" \
                                            % (x, y, self._cursor._series, z_obs))
            else:
                z = (x/self._cursor._xmean.to(sess.spec._xunit).value-1)#*(1+sess.spec._rfz)
                #print(z)
                for c, xem in zip(self._cursor_line, self._cursor._xem):
                    c.set_xdata((xem*(1+z)*au.nm).to(sess.spec._xunit))
                    c.set_alpha(0.5)
                #self._canvas.draw()
                z_obs = z*(1+sess.spec._rfz)
                focus._textbar.SetLabel("x=%2.4f, y=%2.4e; z[%s]=%2.5f" \
                                        % (x, y, self._cursor._series, z_obs))
            self._cursor._z = z
        else:
            focus._textbar.SetLabel("x=%2.4f, y=%2.4e" % (x, y))
        for l in self._cursor_lines:
            for b in l:
                self._ax.draw_artist(b)

        # Copied from https://matplotlib.org/_modules/matplotlib/backends/backend_wxagg.html
        #self._canvas.update()
        self._canvas.draw_idle()
        #self._canvas.flush_events()

    def _refresh(self, sess, logx=False, logy=False, norm=False, legend=None,
                 xlim=None, ylim=None, title=None, text=None, init_cursor=False):
        sess = np.array(sess, ndmin=1)
        #import datetime as dt
        #start = dt.datetime.now()
        self._text = text
        self._ax.clear()
        self._ax.grid(True, which='both', linestyle=':')
        if title != None:
            self._ax.set_title(title)
        if text != None:
            self._ax.text(0.05, 0.1, text, transform=self._ax.transAxes)
        if init_cursor:
            self._cursor_lines = []

        cmc = plt.cm.get_cmap('tab10').colors
        self._canvas_dict = {#'spec_x_y': (GraphSpectrumXY,0,1.0),
                             #'spec_x_y_det': (GraphSpectrumXYDetail,0],1.0),
                             #'spec_x_dy': (GraphSpectrumXDy,0,0.5),
                             #'spec_x_conv': (GraphSpectrumXConv,3,0.5),
                             #'lines_x_y': (GraphLineListXY,2,1.0),
                             #'spec_x_ylinemask': (GraphSpectrumXYLinesMask,2,0.5),
                             #'spec_nodes_x_y': (GraphSpectrumNodesXY,1,1.0),
                             #'spec_x_cont': (GraphSpectrumXCont,8,1.0),
                             #'spec_form_x': (GraphSpectrumFormX,7,0.5),
                             #'spec_x_model': (GraphSpectrumXModel,9,1.0),
                             #'spec_x_yfitmask': (GraphSpectrumXYFitMask,9,0.5),
                             #'spec_x_deabs': (GraphSpectrumXDeabs,9,0.5),
                             #'systs_z_series': (GraphSystListZSeries,2,1.0),
                             'cursor_z_series': (GraphCursorZSeries,3,0.5),
                             'spec_h2o_reg': (GraphSpectrumH2ORegion,4,0.15)
                             }


        #print([self._gui._sess_sel== i for i in self._gui._sess_items])
        c_index = [self._canvas_dict[s][1]\
        #c_index = [(self._canvas_dict[s][1]\
                    #+max(len(self._gui._sess_item_sel),1)-1) % 10 \
                   for s in self._sel]
        self._canvas_l = [self._canvas_dict[s][0] for s in self._sel]
        #if len(self._gui._sess_item_sel) < 2:
        self._color_l = [cmc[i] for i in c_index]
        self._alpha_l = [self._canvas_dict[s][2] for s in self._sel]

        # First selected session sets the units of the axes
        self._xunit = GraphSpectrumXY(sess[0])._x.unit
        self._yunit = GraphSpectrumXY(sess[0], norm)._y.unit
        self._ax.set_xlabel(self._xunit)
        self._ax.set_ylabel(self._yunit)

        #self._ax.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.2e}'))

        #f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
        #g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
        #self._ax.yaxis.set_major_formatter(mticker.FuncFormatter(g))

        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        self._ax.yaxis.set_major_formatter(yfmt)


        # Rest frame axis
        if sess[0].spec._rfz != 0.0:
            self._ax.set_xlabel(str(self._xunit)+", rest frame (z = %3.3f)"
                                % sess[0].spec._rfz)
            if self._axt == None:
                self._axt = self._ax.twiny()
            try:
                self._axt.set_xlabel(str(self._gui._sess_sel.spec._xunit))
            except:
                self._axt.set_xlabel(str(self._xunit))
            self._axt_mode = 'rf'
        else:
            try:
                self._axt.remove()
            except:
                pass
            self._axt = None
        #self._c = 0  # Color

        # Redshift axis
        if hasattr(sess[0], '_ztrans'):
            self._axt = self._ax.twiny()
            self._axt.set_xlabel('%s redshift' % sess[0]._ztrans)
            self._axt_mode = 'z'
        else:
            try:
                self._axt.remove()
            except:
                pass
            self._axt = None

        if logx:
            self._ax.set_xscale('log')
            try:
                self._axt.set_xscale('log')
            except:
                pass
        if logy:
            self._ax.set_yscale('log')
            try:
                self._axt.set_yscale('log')
            except:
                pass

        autoxlim = False
        autoylim = False
        if self == self._gui._graph_main._graph:
            for l in self._gui._graph_main._lim.split('\n'):
                ls = l.split('=')
                if ls[-1][0]=='(' and ls[-1][-1]==')':
                    if ls[0]=='xlim':
                        try:
                            xlim = tuple(map(float, ls[-1][1:-1].split(',')))
                        except:
                            pass
                    if ls[0]=='ylim':
                        try:
                            ylim = tuple(map(float, ls[-1][1:-1].split(',')))
                        except:
                            pass
                elif ls[-1]!='auto':
                    logging.error(msg_lim(ls[0]))
                """
                else:
                    if ls[0]=='xlim':
                        autoxlim = True
                    if ls[0]=='ylim':
                        autoylim = True
                """

        if xlim is not None and not autoxlim:
            self._ax.set_xlim(xlim)
        if ylim is not None and not autoylim:
            self._ax.set_ylim(ylim)

        for s in [sess[0]]:
            self._seq(s, norm, init_cursor=init_cursor)
        if legend:
            self._ax.legend()

        for s in [sess[0]]:
            if s._shade:
                x = self._gui._sess_sel.spec.x.value
                trans = transforms.blended_transform_factory(
                            self._ax.transData, self._ax.transAxes)
                self._ax.fill_between(x, 0, 1, where=s._shade_where,
                                      transform=trans, color='C1', alpha=0.2)


        self._canvas.draw()
        #self._canvas.flush_events()
        #print(dt.datetime.now()-start)

        self._ax.callbacks.connect('xlim_changed', self._on_zoom)
        self._ax.callbacks.connect('ylim_changed', self._on_zoom)
        dl = self._ax.__dict__['dataLim']
        #print(dl)
        self._gui._data_lim = (dl.x0, dl.x1, dl.y0, dl.y1)
        print(self._gui._data_lim)


    def _reg_shade(self):
        sess = self._gui._sess_sel
        x = sess.spec.x.value

        sess._shade_where = np.logical_and(x>sess._clicks[0][0],
                                           x<sess._clicks[1][0])
        """
        trans = transforms.blended_transform_factory(
                    self._ax.transData, self._ax.transAxes)
        self._ax.fill_between(x, 0, 1, where=self._shade_where,
                              transform=trans,
                              color='C1', alpha=0.2)
        """
        sess._shade = True
        #self._refresh(sess, xlim=self._ax.get_xlim(), ylim=self._ax.get_ylim())

        #x = self._gui._sess_sel.spec.x.value
        trans = transforms.blended_transform_factory(
                    self._ax.transData, self._ax.transAxes)

        shade = self._ax.fill_between(x, 0, 1, where=sess._shade_where,
                                      transform=trans, color='C1', alpha=0.2)

        self._canvas.draw()
        shade.remove()

    def _seq(self, sess, norm, init_cursor=True):
        detail = self._panel != self._gui._graph_main._panel
        if detail:
            focus = self._gui._graph_det
        else:
            focus = self._gui._graph_main
        self._seq_core(sess, norm, init_cursor, detail, focus, self._ax)

    def _seq_core(self, sess, norm, init_cursor, detail, focus, ax,
                  cursor_list=['cursor']):
#        ax = self._ax

        xunit_orig = dc(sess.spec.x.unit)
        #print(detail, sess.spec.x.unit)
        #print(self._zem)
        #if detail: sess.cb.x_convert(zem=self._zem)
        #print(detail, sess.spec.x.unit)
        self._systs_id = False


        for e in focus._elem.split('\n'):
        #for e in self._gui._graph_main._elem.split('\n'):
            try:
                sel, struct, xcol, ycol, mcol, mode, style, width, color, alpha\
                    = e.split(',')
                #print(sel, struct, xcol, ycol, mcol, mode, style, width, color, alpha)
                sessw = self._gui._sess_item_list.index(int(sel))
                sess = self._gui._sess_list[sessw]
                #sess = self._gui._sess_sel
                xunit = sess.spec.x.unit
                if mcol != 'None':
                    label = '%s, %s (%s)' % (struct, ycol, mcol)
                else:
                    label = '%s, %s' % (struct, ycol)

                if struct == 'systs': self._systs_id = True
                if struct in ['spec','lines','nodes','systs','feats']:
                    t = getattr(sess, struct).t
                    if mode != 'axhline':
                        x = dc(t[xcol])
                    if mode != 'axvline':
                        try:
                            y = int(ycol)
                        except:
                            y = dc(t[ycol])
                    if mcol not in ['None', 'none', None]:
                        x[t[mcol]==0] = np.nan
                    if norm and 'cont' in t.colnames and t[ycol].unit == t[ycol].unit and len(y)==len(t['cont']):
                        y = y/t['cont']
                #if detail: print(struct, xcol, x[0])
                #print(sel, struct, xcol, ycol, mcol, mode, style, width, color, alpha)
                if struct in ['systs', 'cursor']:
                    if xcol == 'z' :
                    #if struct == 'systs':
                        z = sess.systs.z
                        series = sess.systs.series
                        z_list = [[zf]*len(trans_parse(s)) for zf,s in zip(z,series)]
                        series_list = [trans_parse(s) for s in series]
                        z_flat = np.array([z for zl in z_list for z in zl])
                        series_flat = np.array([s for sl in series_list for s in sl])
                    else:
                        z = float(xcol)
                        series = sess._cursors[xcol]._series
                        z_flat = np.array([z]*len(trans_parse(series)))
                        series_flat = trans_parse(series)
                    xem = np.array([xem_d[sf].to(au.nm).value \
                                    for sf in series_flat]) * au.nm


                    if hasattr(self._gui._sess_sel.spec, '_rfz'):
                        x = xem*(1+z_flat/(1+self._gui._sess_sel.spec._rfz))
                    else:
                        x = xem*(1+z_flat)

                    # Convert to spectrum x units if needed
                    try:
                        x = x_convert(x, self._gui._sess_sel.spec._zem,
                                      self._gui._sess_sel.spec._xunit)
                        x_iswave = True
                    except:
                        x_iswave = False
                    """
                    try:
                        x = x.to(sess.spec._xunit)
                        #sess.spec_x_convert(self, zem=0, xunit=au.km/au.s)
                        x_iswave = True
                    except:
                        x_iswave = False
                    """
                    #print(graph._xs)
                    #print(self._zems, self._series, self._axes, self._ax)
                    #print(zems)
                    self._systs_l = x
                    if detail:
                        #print(z, self._zem)
                        #print((1+self._zem)*121.567)
                        #print(self._zem)
                        #x = np.log(x.value/((1+self._zem)*121.567))*aconst.c.to(au.km/au.s)
                        for k in self._axes:
                            if self._axes[k] == ax:
                                zem = self._zems[k]
                        if x_iswave:
                            x = np.log(x.to(au.nm).value/((1+zem)*121.567))*aconst.c.to(au.km/au.s)
                        else:
                            x = np.log(x.value/((1+zem)*121.567))*aconst.c.to(au.km/au.s)

                        #print(set(zip(series_flat,x)))
                    self._systs_series = series_flat
                    self._systs_z = z_flat
                    self._systs_x = x

                    if hasattr(self._gui._graph_main, '_z_sel'):
                        z_sel = self._gui._graph_main._z_sel
                        series_sel = self._gui._graph_main._series_sel
                        xem_sel = np.array([xem_d[s].to(au.nm).value \
                                            for s in series_sel])
                        x_sel = xem_sel*(1+z_sel)
                    else:
                        x_sel = []
                try:
                    kwargs = {}
                    if mode in ['plot', 'step', 'axvline']:
                        kwargs['linestyle'] = style
                        kwargs['linewidth'] = width
                    if mode == 'step':
                        kwargs['where'] = 'mid'
                    if mode == 'scatter':
                        kwargs['marker'] = style
                        kwargs['s'] = (5*float(width))**2
                    kwargs['color'] = color
                    kwargs['alpha'] = float(alpha)
                    if mode == 'axvline':
                        #print(self._ax, x.value)
                        for xi in x.value:
                            if xi==x[0].value:
                                getattr(ax, mode)(xi, label='systs', **kwargs)
                            else:
                                getattr(ax, mode)(xi, **kwargs)
                        if focus == self._gui._graph_main:
                            for xi_sel in x_sel:
                                kwargs['linestyle'] = '-'
                                kwargs['linewidth'] = 10.0
                                kwargs['color'] = 'yellow'
                                kwargs['alpha'] = 0.3
                                if xi_sel==x_sel[0]:
                                    getattr(ax, mode)(xi_sel,
                                        label='systs (selected)', **kwargs)
                                else:
                                    getattr(ax, mode)(xi_sel, **kwargs)
                    else:
                        if type(y) in [int, float]:
                            y = [y]*len(x)
                        getattr(ax, mode)(x, y, label=label, **kwargs)

                    if struct in cursor_list:
                        trans = transforms.blended_transform_factory(
                                ax.transData, ax.transAxes)
                        for xi, s, z in zip(x.value, series_flat, z_flat):
                            if xi > ax.get_xlim()[0] \
                                and xi < ax.get_xlim()[1]:
                                kwargs_text = {}
                                kwargs_text['color'] = color
                                kwargs_text['alpha'] = float(alpha)
                                kwargs_text['size'] = (float(width)+1)*5
                                kwargs_text['transform'] = trans
                                kwargs_text['rotation'] = 90
                                kwargs_text['ha'] = 'right'
                                kwargs_text['va'] = 'bottom'
                                if hasattr(self._gui._sess_sel.spec, '_rfz_man'):
                                    zz = self._gui._sess_sel.spec._rfz_man
                                    z = z*(1+zz)+zz
                                elif hasattr(self._gui._sess_sel.spec, '_rfz'):
                                    z += self._gui._sess_sel.spec._rfz

                                if z > 1e-10:
                                    ax.text(xi, 0.05, s, **kwargs_text)
                                    kwargs_text['va'] = 'top'
                                    ax.text(xi, 0.95, "%3.5f" % z, **kwargs_text)
                                else:
                                    kwargs_text['va'] = 'top'
                                    ax.text(xi, 0.95, s, **kwargs_text)
                    if struct == 'systs':
                        self._systs_color = color
                except:
                    logging.error("I can't parse this graph specification: %s." % e)

            except:
                pass

        if hasattr(sess.spec, '_stats_text_red'):
            ax.text(0.98, 0.95, sess.spec._stats_text_red, va='top', ha='right',
                          transform=ax.transAxes)

        self._check_units(sess, 'x')
        self._check_units(sess, 'y')
        for z, (s, c, a) \
            in enumerate(zip(self._canvas_l, self._color_l, self._alpha_l)):
            try:
                gs = s(sess, norm)
                if gs._type == 'axvline':
                    for x in gs._x:
                        #print(x)
                        #print(**gs._kwargs)
                        ax.axvline(x.to(self._xunit).value,
                                         color=c, alpha=a, linewidth=1.5,
                                         **gs._kwargs)
                        gs._kwargs.pop('label', None)
                    try:
                        for x in gs._xalt:
                            #print(x)
                            ax.axvline(x.to(self._xunit).value,
                                             color=c, alpha=a,
                                             linewidth=0.5, **gs._kwargs)
                            gs._kwargs.pop('label', None)

                    except:
                        pass
                elif gs._type == 'fill_between':
                    ylim = ax.get_ylim()
                    trans = transforms.blended_transform_factory(
                                ax.transData, ax.transAxes)
                    if hasattr(sess.spec, '_rfz_man'):
                        reg = h2o_reg/(1+sess.spec._rfz_man)
                    else:
                        reg = h2o_reg/(1+sess.spec._rfz)
                    x = gs._x.to(au.nm).value
                    where = np.logical_and(x>reg[0][0], x<reg[0][1])\
                                +np.logical_and(x>reg[1][0], x<reg[1][1])\
                                +np.logical_and(x>reg[2][0], x<reg[2][1])
                    ax.fill_between(gs._x.to(self._xunit).value, 0, 1,
                                          where=where,
                                          transform=trans, color='gray', alpha=a)
                    ax.set_ylim(ylim)

                elif gs._type == 'text':
                    for (x, t) in zip(gs._x, gs._y):
                        #print(x,t)
                        ax.text(x.to(self._xunit).value, 0.8, t,
                                      horizontalalignment='center',
                                      transform=trans)
                        #print("here", x,t)
                        gs._kwargs.pop('label', None)
                        #print("after", x,t)
                elif gs._type == 'axvline_special':
                    try:
                        cz = self._cursor._z
                    except:
                        pass

                    self._cursor = gs

                    self._cursor_line = []
                    if gs._z == 0:
                        #break
                        #gs._z = (1+self._zems[self._text])*xem_d['Ly_a']/gs._xmean-1
                        #print('gs_xem', gs._xem)
                        #gs._z = (1+self._zem)*xem_d['Ly_a']/(np.min(gs._xem)*au.nm)-1
                        try:
                            gs._z = cz
                        except:
                            gs._z = self._z
                        gs._x = gs._xem*(1+gs._z)*au.nm
                        xem = self._xs[self._text]
                        #xem = self._x
                        #print(xem)
                        #equiv = [(au.nm, au.km/au.s,
                        #         lambda x: np.log(x/xem.value)*aconst.c.to(au.km/au.s),
                        #         lambda x: np.exp(x/aconst.c.to(au.km/au.s).value)*xem.value)]
                        gs._x = np.log(gs._x/xem)*aconst.c.to(au.km/au.s)
                    for i, x in enumerate(gs._x):
                        if i==1:
                            del gs._kwargs['label']
                        self._cursor_line.append(
                            ax.axvline(
                                x.to(self._xunit).value, #alpha=0,
                                color=c, alpha=a, linewidth=2,
                                **gs._kwargs))
                        """
                        if focus==self._gui._graph_main:
                            ax.axvline(
                                x.to(self._xunit).value,
                                color='C3', alpha=0.3, linewidth=10)
                        """

                    self._cursor_lines.append(self._cursor_line)
                else:
                    graph = getattr(ax, gs._type)
                    graph(gs._x, gs._y, zorder=z, color=c, alpha=a,
                          **gs._kwargs)
                #self._c += 1
            except:
                pass
            if self._axt != None:
                if self._axt_mode == 'rf':
                    self._axt.set_xlim(np.array(ax.get_xlim()) \
                                       * (1+sess.spec._rfz))
                if self._axt_mode == 'z':
                    self._axt.set_xlim((np.array(ax.get_xlim())*self._xunit \
                        / xem_d[sess._ztrans]).to(au.dimensionless_unscaled)-1)

        #if detail: sess.cb.x_convert(zem=self._zem, xunit=xunit_orig)

        """
        for c in self._gui._graph_main._cols_sel.split(','):
            if c in sess.spec.t.colnames:
                y = sess.spec.t[c]
                if norm and 'cont' in sess.spec._t.colnames:
                    y = y/sess.spec.t['cont']
                ax.plot(sess.spec.x, y)
        """


class GraphLineListXY(object):
    def __init__(self, sess, norm=False):
        self._type = 'scatter'
        self._x = sess.lines.x
        self._y = sess.lines.y
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.lines._t['cont']
        self._kwargs = {'marker':'+', 'label':sess.name+", lines"}


class GraphSpectrumFormX(object):
    def __init__(self, sess):
        self._type = 'axvline'
        self._x = sess.spec_form.x
        self._kwargs = {'linewidth':1.0, 'linestyle': ':',
                        'label':sess.spec._meta['instr']+" spectral format"}

class GraphSpectrumNodesXY(object):
    def __init__(self, sess, norm=False):
        self._type = 'scatter'
        self._x = sess.nodes.x
        self._y = sess.nodes.y
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.nodes._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", nodes"}

class GraphSpectrumXY(object):
    def __init__(self, sess, norm=False):
        self._type = 'step'
        #self._type = 'plot'
        self._x = sess.spec.x
        self._y = sess.spec.y
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name, 'where':'mid'}
        #self._kwargs = {'lw':1.0, 'label':sess.name}


class GraphSpectrumH2ORegion(GraphSpectrumXY):
    def __init__(self, sess, norm=False):
        super(GraphSpectrumH2ORegion, self).__init__(sess)
        self._type = 'fill_between'


class GraphSpectrumXCont(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXCont, self).__init__(sess)
        self._type = 'plot'
        self._y = sess.spec._t['cont']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", continuum"}

class GraphSpectrumXConv(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXConv, self).__init__(sess)
        self._y = sess.spec._t['y_conv']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", convolved"}

class GraphSpectrumXDy(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXDy, self).__init__(sess)
        self._y = sess.spec.dy
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", error", 'where':'mid'}

class GraphSpectrumXModel(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXModel, self).__init__(sess)
        self._type = 'plot'
        self._x = sess.spec.x
        self._y = sess.spec._t['model']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", model"}

class GraphSpectrumXDeabs(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXDeabs, self).__init__(sess)
        self._y = sess.spec._t['deabs']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", de-absorbed",
                        'where':'mid'}

class GraphSpectrumXYDetail(object):
    def __init__(self, sess, norm=False):
        self._type = 'scatter'
        self._x = sess._xdet
        self._y = sess._ydet
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'marker':'o'}

class GraphSpectrumXYFitMask(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXYFitMask, self).__init__(sess)
        self._type = 'plot'
        self._x[sess.spec._t['fit_mask']==0] = np.nan
        self._y = sess.spec._t['model']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':3.0, 'label':sess.name+", masked for fitting"}

class GraphSpectrumXYLinesMask(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXYLinesMask, self).__init__(sess)
        self._type = 'step'
        self._x[sess.spec._t['lines_mask']] = np.nan
        try:
            self._y = sess.spec._t['deabs']
        except:
            self._y = sess.spec.y
        self._kwargs = {'label':sess.name+", masked for lines", 'where':'mid'}

class GraphSystListZSeries(object):
    def __init__(self, sess, norm=False):
        #self._type = 'text'
        self._type = 'axvline'
        #"""
        z = sess.systs.z
        series = sess.systs.series

        #z_list = [[zf]*len(series_d[s]) for zf,s in zip(z,series)]
        z_list = [[zf]*len(trans_parse(s)) for zf,s in zip(z,series)]

        #z_flat = np.ravel([[zf]*len(series_d[s]) for zf,s in zip(z,series)])
        #series_list = [series_d[s] for s in series]
        series_list = [trans_parse(s) for s in series]

        #series_flat = np.ravel([series_d[s] for s in series])

        z_flat = np.array([z for zl in z_list for z in zl])
        series_flat = np.array([s for sl in series_list for s in sl])

        kn = np.where(series_flat != 'unknown')
        unkn = np.where(series_flat == 'unknown')

        z_kn = z_flat[kn]
        series_kn = series_flat[kn]
        z_unkn = z_flat[unkn]
        series_unkn = series_flat[unkn]
        #print(z_flat[unkn])
        #print(series_flat[unkn])
        xem_kn = np.array([xem_d[sf].to(au.nm).value for sf in series_kn])\
                   *au.nm
        #print(xem_kn)
        self._x = (1.+z_kn)*xem_kn
        self._xalt = z_unkn * au.nm
        #print(self._x)
        #print(self._xalt)
        self._kwargs = {'linestyle': ':',
                        'label':sess.name+", system components"}
        """
        self._y = series_flat
        series = np.array([sess.systs.series[i[0]]
                           for i in sess.systs._mods_t['id']])
        n = np.array([len(i) for i in sess.systs._mods_t['id']])
        z = np.array(sess.systs._mods_t['z0'])
        xem = np.array([xem_d[s].to(au.nm).value for s in series])*au.nm
        self._x = (1.+z)*xem
        print(self._x)
        self._y = np.array(["%s (%i)\n%3.5f" % (si,ni,zi)
                           for si,ni,zi in zip(series,n,z)])

        self._kwargs = {'marker':'+', 'label':sess.name+", systs"}
        """

class GraphCursorZSeries(object):
    def __init__(self, sess, norm=False):
        self._type = 'axvline_special'
        if hasattr(sess, '_series_sel'):
            self._series = sess._series_sel
        else:
            self._series = 'CIV'
        #self._xem = np.array([xem_d[s].to(au.nm).value \
        #                      for s in series_d[self._series]])
        self._xem = np.array([xem_d[t].value for t in trans_parse(self._series)])
        self._xmean = np.mean(self._xem)*au.nm
        try:
            self._z = np.mean(sess.spec.x).to(au.nm).value/self._xmean.value-1.0
            self._z = sess._z_sel
        except:
            self._z = 0
        """
        if hasattr(sess, '_hwin_sel'):
            self._hwin = sess._hwin_sel
        else:
            self._hwin = 250.0
        """
        self._x = self._xem*(1+self._z)*au.nm
        self._kwargs = {'label':self._series, 'linestyle': '--'}
