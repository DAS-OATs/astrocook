from .frame import Frame
from .vars import *
from astropy import units as au
import numpy as np

class LineList(Frame):
    """Class for line lists

    A LineList is a Frame with methods for handling spectral lines."""

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

    def _copy(self, sel=None):
        copy = super(LineList, self)._copy(sel)
        cols = [c for c in self._t.colnames \
                if c not in ['x', 'xmin', 'xmax', 'y', 'dy']]
        for c in cols:
            copy._t[c] = self._t[c][sel]
        return copy

    """
    def _ew(self, spec, thres=1e-1):

        for l in self._t:
            xmin = float(round(l['xmin'],3))
            xmax = float(round(l['xmax'],3))
            c = np.where(np.spec.l['xmin'])
            spec_reg = self.acs.spec.extract_reg(xmin, xmax)
            cont_reg = self.acs.cont.extract_reg(xmin, xmax)
            l['ew'] = np.sum((1-spec_reg.t['Y']/cont_reg.t['Y'])\
                        *(spec_reg.t['XMAX']-spec_reg.t['XMIN']))\
            *spec_reg.t['X'].unit
            l['EW'] = ew.value
            if ew.value > thres:
                l['XMIN'] = l['X']-0.5*ew.value
                l['XMAX'] = l['X']+0.5*ew.value
    """

    def _syst_cand(self, series, z_start, z_end, dz):

        # Compute all possible redshifts
        z_all = np.ravel([[(x.to(au.nm)/xem_d[t].to(au.nm)).value-1. \
                            for x in self.x] for t in series_d[series]])

        # Select values within [z_start, z_end]
        (z_min, z_max) = (z_start, z_end) if z_start < z_end \
            else (z_end, z_start)
        z_sel = z_all[np.logical_and(z_all>z_min, z_all<z_max)]

        # Find coincidences
        if len(series_d[series]) > 1:
            z_sort = np.sort(np.ravel(z_sel))
            z_range = z_sort[np.where(np.ediff1d(z_sort)<dz)]
            z_range = z_range if z_start < z_end else z_range[::-1]
        else:
            z_range = z_sel

        return z_range
