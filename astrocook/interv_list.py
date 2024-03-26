from astropy import units as au
from astropy import table as at


class Interv(object):
    """ Class for intervals """

    def __init__(self, xmin, xmax, xunit=au.nm):
        self._xmin = xmin
        self._xmax = xmax
        self._xunit = xunit


class IntervList(object):
    """ Class for interval lists """


    def __init__(self, xmin=[], xmax=[], xunit=au.nm):
        self._xmin = xmin
        self._xmax = xmax
        self._xunit = xunit

        t = at.Table()
        t['xmin'] = at.Column(np.array(xmin, ndmin=1), dtype=dtype, unit=xunit)
        t['xmax'] = at.Column(np.array(xmax, ndmin=1), dtype=dtype, unit=xunit)
        self._t = t
        self._l = []

    @property
    def t(self):
        return self._t


    def _add(self, xmin, xmax, xunit=au.nm):
        interv = Interv(xmin, xmax, xunit)
        self._t.add_row([xmin, xmax])
        self._l.append(interv)
        return 0
