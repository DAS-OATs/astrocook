from .frame import Frame
from astropy import units as au

class LineList(Frame):
    """Class for line lists

    A Line List is a Frame with methods for handling spectral lines."""

    def __init__(self,
                 x=[],
                 xmin=[],
                 xmax=[],
                 y=[],
                 dy=[],
                 xunit=au.nm,
                 yunit=au.erg/au.cm**2/au.s/au.nm,
                 meta={},
                 dtype=float):
        super(LineList, self).__init__(x, xmin, xmax, y, dy, xunit, yunit, meta,
                                       dtype)
