from .message import *
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
import numpy as np

prefix = "Graph:"

class Graph(object):

    def __init__(self, panel, gui, sel):
        self._gui = gui
        self._sel = sel
        self._fig = Figure()
        self._ax = self._fig.add_subplot(111)
        #self._fig.tight_layout()#rect=[-0.03, 0.02, 1.03, 1])
        self._plot = FigureCanvasWxAgg(panel, -1, self._fig)
        self._toolbar = NavigationToolbar2WxAgg(self._plot)
        self._toolbar.Realize()

    def _check_units(self, sess, axis='x'):
        unit = axis+'unit'
        _unit = '_'+unit
        if getattr(getattr(sess, 'spec'), _unit) != getattr(self, _unit):
            print(prefix, "I'm converting the %s unit of %s to plot it over the "
                  "data already present." % (axis, sess.name))
            getattr(sess, 'convert_'+axis)(**{unit: getattr(self, _unit)})
            self._gui._panel_sess._refresh()

    def _refresh(self, sess, logx=False, logy=False):
        sess = np.array(sess, ndmin=1)
        self._ax.clear()
        self._plot_dict = {'spec_x_y': GraphSpectrumXY,
                           'spec_x_dy': GraphSpectrumXDy,
                           'spec_x_conv': GraphSpectrumXConv,
                           'lines_x_y': GraphLineListXY,
                           'spec_x_ymask': GraphSpectrumXYMask,
                           'spec_nodes_x_y': GraphSpectrumNodesXY,
                           'spec_x_cont': GraphSpectrumXCont,
                           'spec_form_x': GraphSpectrumFormX}
        self._plot_list = [self._plot_dict[s] for s in self._sel]

        # First selected session sets the units of the axes
        self._xunit = GraphSpectrumXY(sess[0])._x.unit
        self._yunit = GraphSpectrumXY(sess[0])._y.unit
        self._ax.set_xlabel(self._xunit)
        self._ax.set_ylabel(self._yunit)
        self._c = 0  # Color
        if logx:
            self._ax.set_xscale('log')
        if logy:
            self._ax.set_yscale('log')

        for s in sess:
            self._seq(s)
        self._ax.legend()
        self._plot.draw()

    def _seq(self, sess):
        self._check_units(sess, 'x')
        self._check_units(sess, 'y')
        for z, s in enumerate(self._plot_list):
            try:
                gs = s(sess)
                if gs._type == 'axvline':
                    for x in gs._x:
                        self._ax.axvline(x.to(self._xunit).value,
                                         color='C'+str(self._c), **gs._kwargs)
                        gs._kwargs.pop('label', None)
                else:
                    graph = getattr(self._ax, gs._type)
                    graph(gs._x, gs._y, zorder=z, color='C'+str(self._c),
                          **gs._kwargs)
                self._c += 1
            except:
                pass

class GraphLineListXY(object):
    def __init__(self, sess):
        self._type = 'scatter'
        self._x = sess.lines.x
        self._y = sess.lines.y
        self._kwargs = {'marker':'+', 'label':sess.name+", lines"}

class GraphSpectrumFormX(object):
    def __init__(self, sess):
        self._type = 'axvline'
        self._x = sess.spec_form.x
        self._kwargs = {'lw':1.0, 'linestyle': ':',
                        'label':sess.spec._meta['instr']+" spectral format"}

class GraphSpectrumNodesXY(object):
    def __init__(self, sess):
        self._type = 'plot'
        self._x = sess.nodes.x
        self._y = sess.nodes.y
        self._kwargs = {'lw':1.0, 'label':sess.name+", nodes"}

class GraphSpectrumXY(object):
    def __init__(self, sess):
        self._type = 'plot'
        self._x = sess.spec.x
        self._y = sess.spec.y
        self._kwargs = {'lw':1.0, 'label':sess.name}

class GraphSpectrumXCont(GraphSpectrumXY):

    def __init__(self, sess):
        super(GraphSpectrumXCont, self).__init__(sess)
        self._y = sess.spec._t['cont']
        self._kwargs = {'lw':1.0, 'label':sess.name+", continuum"}

class GraphSpectrumXConv(GraphSpectrumXY):

    def __init__(self, sess):
        super(GraphSpectrumXConv, self).__init__(sess)
        self._y = sess.spec._t['conv']
        self._kwargs = {'lw':1.0, 'label':sess.name+", convolved"}

class GraphSpectrumXDy(GraphSpectrumXY):

    def __init__(self, sess):
        super(GraphSpectrumXDy, self).__init__(sess)
        self._y = sess.spec.dy
        self._kwargs = {'lw':1.0, 'label':sess.name+", error"}

class GraphSpectrumXYMask(GraphSpectrumXY):

    def __init__(self, sess):
        super(GraphSpectrumXYMask, self).__init__(sess)
        self._x[sess.spec._t['lines_mask']] = np.nan
        self._kwargs = {'lw':1.0, 'label':sess.name+", masked"}
