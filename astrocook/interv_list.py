from .message import msg_param_swap
from astropy import units as au
from astropy import table as at
import numpy as np

class Interv(object):
    """ Class for intervals """

    def __init__(self, xmin, xmax, xunit=au.nm):
        self._xmin = xmin
        self._xmax = xmax
        self._xunit = xunit


class IntervList(object):
    """ Class for interval lists """


    def __init__(self, xmin=[], xmax=[], xunit=au.nm, dtype=float):
        self._xmin = xmin
        self._xmax = xmax
        self._xunit = xunit

        t = at.Table()
        t['xmin'] = at.Column(np.array(xmin, ndmin=1), dtype=dtype, unit=xunit)
        t['xmax'] = at.Column(np.array(xmax, ndmin=1), dtype=dtype, unit=xunit)
        self._t = t
        self._l = []
        for m, M in zip(xmin, xmax):
            interv = Interv(xmin, xmax, xunit)
            self._l.append(interv)


    @property
    def t(self):
        return self._t


    def _add(self, xmin, xmax, xunit=au.nm):
        if xmin > xmax:
            xmin, xmax = xmax, xmin
            logging.warning(msg_param_swap)

        interv = Interv(xmin, xmax, xunit)
        self._t.add_row([xmin, xmax])
        self._l.append(interv)
        return 0


    def _clear(self):
        self.__init__()
