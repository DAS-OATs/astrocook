from .message import *
from .vars import *
from astropy import units as au
#from copy import deepcopy as dc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
from matplotlib.figure import Figure
import matplotlib.transforms as transforms
import numpy as np

prefix = "Graph:"

class Graph(object):

    def __init__(self, panel, gui, sel):
        self._gui = gui
        self._sel = sel
        self._fig = Figure()
        self._ax = self._fig.add_subplot(111)
        self._c = 0
        #self._fig.tight_layout()#rect=[-0.03, 0.02, 1.03, 1])
        self._plot = FigureCanvasWxAgg(panel, -1, self._fig)
        self._toolbar = NavigationToolbar2WxAgg(self._plot)
        self._toolbar.Realize()

    def _check_units(self, sess, axis='x'):
        unit = axis+'unit'
        _unit = '_'+unit
        if getattr(getattr(sess, 'spec'), _unit) != getattr(self, _unit) \
            and getattr(self, _unit) != au.dimensionless_unscaled:
            print(prefix, "I'm converting the %s unit of %s to plot it over the "
                  "data already present." % (axis, sess.name))
            getattr(sess, 'convert_'+axis)(**{unit: getattr(self, _unit)})
            self._gui._panel_sess._refresh()

    def _refresh(self, sess, logx=False, logy=False, norm=False):
        sess = np.array(sess, ndmin=1)
        self._ax.clear()
        self._plot_dict = {'spec_x_y': GraphSpectrumXY,
                           'spec_x_dy': GraphSpectrumXDy,
                           'spec_x_conv': GraphSpectrumXConv,
                           'lines_x_y': GraphLineListXY,
                           'spec_x_ymask': GraphSpectrumXYMask,
                           'spec_nodes_x_y': GraphSpectrumNodesXY,
                           'spec_x_cont': GraphSpectrumXCont,
                           'spec_form_x': GraphSpectrumFormX,
                           'spec_x_model': GraphSpectrumXModel,
                           'spec_x_deabs': GraphSpectrumXDeabs,
                           'systs_z_series': GraphSystListZSeries}
        #print(self._sel)
        self._plot_list = [self._plot_dict[s] for s in self._sel]

        # First selected session sets the units of the axes
        self._xunit = GraphSpectrumXY(sess[0])._x.unit
        self._yunit = GraphSpectrumXY(sess[0], norm)._y.unit
        self._ax.set_xlabel(self._xunit)
        self._ax.set_ylabel(self._yunit)
        if sess[0].spec._rfz != 0.0:
            self._ax.set_xlabel(str(self._xunit)+", rest frame (z = %3.2f)"
                                % sess[0].spec._rfz)
        self._c = 0  # Color
        if logx:
            self._ax.set_xscale('log')
        if logy:
            self._ax.set_yscale('log')

        for s in sess:
            self._seq(s, norm)
        self._ax.legend()
        self._plot.draw()

    def _seq(self, sess, norm):
        self._check_units(sess, 'x')
        self._check_units(sess, 'y')
        for z, s in enumerate(self._plot_list):
            try:
                gs = s(sess, norm)
                #print(gs._type)
                if gs._type == 'axvline':
                    for x in gs._x:
                        #print(x)
                        #print(**gs._kwargs)
                        self._ax.axvline(x.to(self._xunit).value,
                                         color='C'+str(self._c), linewidth=1.5,
                                         **gs._kwargs)
                        gs._kwargs.pop('label', None)
                    try:
                        for x in gs._xalt:
                            #print(x)
                            self._ax.axvline(x.to(self._xunit).value,
                                             color='C'+str(self._c),
                                             linewidth=0.5, **gs._kwargs)
                            gs._kwargs.pop('label', None)

                    except:
                        pass
                elif gs._type == 'text':
                    trans = transforms.blended_transform_factory(
                                self._ax.transData, self._ax.transAxes)
                    for (x, t) in zip(gs._x, gs._y):
                        #print(x,t)
                        self._ax.text(x.to(self._xunit).value, 0.8, t,
                                      horizontalalignment='center',
                                      transform=trans)
                        #print("here", x,t)
                        gs._kwargs.pop('label', None)
                        #print("after", x,t)
                else:
                    graph = getattr(self._ax, gs._type)
                    graph(gs._x, gs._y, zorder=z, color='C'+str(self._c),
                          **gs._kwargs)
                self._c += 1
            except:
                pass

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
        self._type = 'plot'
        self._x = sess.nodes.x
        self._y = sess.nodes.y
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.nodes._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", nodes"}

class GraphSpectrumXY(object):
    def __init__(self, sess, norm=False):
        self._type = 'plot'
        self._x = sess.spec.x
        self._y = sess.spec.y
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name}

class GraphSpectrumXCont(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXCont, self).__init__(sess)
        self._y = sess.spec._t['cont']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", continuum"}

class GraphSpectrumXConv(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXConv, self).__init__(sess)
        self._y = sess.spec._t['conv']
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", convolved"}

class GraphSpectrumXDy(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXDy, self).__init__(sess)
        self._y = sess.spec.dy
        if norm and 'cont' in sess.spec._t.colnames:
            self._y = self._y/sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", error"}

class GraphSpectrumXModel(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXModel, self).__init__(sess)
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
        self._kwargs = {'lw':1.0, 'label':sess.name+", de-absorbed"}

class GraphSpectrumXYMask(GraphSpectrumXY):

    def __init__(self, sess, norm=False):
        super(GraphSpectrumXYMask, self).__init__(sess)
        self._type = 'scatter'
        self._x[sess.spec._t['lines_mask']] = np.nan
        self._kwargs = {'label':sess.name+", masked"}

class GraphSystListZSeries(object):
    def __init__(self, sess, norm=False):
        #self._type = 'text'
        self._type = 'axvline'
        #"""
        z = sess.systs.z
        series = sess.systs.series

        z_list = [[zf]*len(series_d[s]) for zf,s in zip(z,series)]

        #z_flat = np.ravel([[zf]*len(series_d[s]) for zf,s in zip(z,series)])
        series_list = [series_d[s] for s in series]
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
                        'label':sess.name+", systs components"}
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
