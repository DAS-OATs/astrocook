from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg

seq = [('plot', 'spec', 'x', 'y', 'C0', 'Spectrum.y'),
       ('plot', 'spec', 'x', 'conv', 'C1', 'Spectrum.conv'),
       ('scatter', 'lines', 'x', 'y', 'C2', 'Line.y')]


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
        tab = getattr(sess, seq[0][1]).t
        self._ax.set_xlabel(tab[seq[0][2]].unit)
        self._ax.set_ylabel(tab[seq[0][3]].unit)
        for s in seq:
            #tab = None
            if hasattr(sess, s[1]):
                try:
                    tab = getattr(sess, s[1]).t
                    mode = getattr(self._ax, s[0])
                    mode(tab[s[2]], tab[s[3]], color=s[4], label=s[5])
                except:
                    pass
