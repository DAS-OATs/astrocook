from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg
import numpy as np

class Graph(object):

    def __init__(self, panel):
        self._fig = Figure()
        self._ax = self._fig.add_subplot(111)
        #self._fig.tight_layout(rect=[-0.03, 0.02, 1.03, 1])
        self._plot = FigureCanvasWxAgg(panel, -1, self._fig)
        self._toolbar = NavigationToolbar2WxAgg(self._plot)
        self._toolbar.Realize()

    def _refresh(self, sess):
        self._ax.clear()
        self._seq(sess)
        self._ax.legend()
        self._plot.draw()

    def _seq(self, sess):
        self._ax.set_xlabel(GraphSpectrum(sess)._x.unit)
        self._ax.set_ylabel(GraphSpectrum(sess)._y.unit)
        for z, s in enumerate(seq):
            try:
                gs = s(sess)
                graph = getattr(self._ax, gs._type)
                graph(gs._x, gs._y, zorder=z, **gs._kwargs)
            except:
                pass

class GraphLineList(object):
    def __init__(self, sess):
        self._type = 'scatter'
        self._x = sess.lines.x
        self._y = sess.lines.y
        self._kwargs = {'color': 'C2', 'label': sess.name+", lines"}

class GraphSpectrum(object):

    def __init__(self, sess):
        self._type = 'plot'
        self._x = sess.spec.x
        self._y = sess.spec.y
        self._kwargs = {'color': 'C0', 'label': sess.name}

class GraphSpectrumCont(GraphSpectrum):

    def __init__(self, sess):
        super(GraphSpectrumCont, self).__init__(sess)
        self._y = sess.spec._t['cont']
        self._kwargs = {'color': 'C4', 'label': sess.name+", continuum"}

class GraphSpectrumConv(GraphSpectrum):

    def __init__(self, sess):
        super(GraphSpectrumConv, self).__init__(sess)
        self._y = sess.spec._t['conv']
        self._kwargs = {'color': 'C1', 'label': sess.name+", convolved"}

class GraphSpectrumMask(GraphSpectrum):

    def __init__(self, sess):
        super(GraphSpectrumMask, self).__init__(sess)
        self._x[sess.spec._t['lines_mask']] = np.nan
        self._kwargs = {'color': 'C3', 'label': sess.name+", masked"}

seq = [GraphSpectrum, GraphSpectrumConv, GraphLineList, GraphSpectrumMask,
       GraphSpectrumCont]
