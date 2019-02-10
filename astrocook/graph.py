from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, \
    NavigationToolbar2WxAgg

seq = [('plot', 'spec', 'x', 'y', 'C0', 'y'),
       ('plot', 'spec', 'x', 'conv', 'C1', 'conv'),
       ('scatter', ('spec', '_peaks'), 'x', 'conv', 'C2', 'peaks')]


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
        self._try(sess)
        self._ax.legend()
        self._plot.draw()

    def _try(self, sess):
        tab = getattr(sess, seq[0][1]).t
        self._ax.set_xlabel(tab[seq[0][2]].unit)
        self._ax.set_ylabel(tab[seq[0][3]].unit)
        for s in seq:
            mode = getattr(self._ax, s[0])
            try:
                if len(s[1]) == 1:
                    tab = getattr(sess, s[1]).t
                elif len(s[1]) == 2:
                    tab = getattr(getattr(sess, s[1][0]), s[1][1]).t
                mode(tab[s[2]], tab[s[3]], color=s[4], label=s[5])
            except:
                pass
