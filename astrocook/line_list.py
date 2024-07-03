from .frame import Frame
from .functions import *
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
                 cont=None,
                 source=None,
                 xunit=au.nm,
                 yunit=au.erg/au.cm**2/au.s/au.nm,
                 meta={},
                 dtype=float,
                 row=None,
                 kind=None):
        if row is not None:
            super(LineList, self).__init__(
                row['x'], row['xmin'], row['xmax'], row['y'], row['dy'],
                xunit, yunit, meta, dtype)
        else:
            super(LineList, self).__init__(
                x, xmin, xmax, y, dy, xunit, yunit, meta, dtype)

        if cont is not None:
            self._t['cont'] = cont
        if source is not None:
            self._t['source'] = source
        if kind is not None:
            self._t['kind'] = kind


    def _add(self, x, xmin, xmax, y, dy, source=None):
        """ @brief Add a system to a system list.
        """

        self._t.add_row([x, xmin, xmax, y, dy, source])

        return 0


    def _clean(self):
        i = 0
        while i < len(self.t)-2:
            l = self._t[i]
            ln = self._t[i+1]
            if ln['x']-0.2*(ln['x']-ln['xmin'])<l['x']:
                l['xmax'] = max(l['xmax'], ln['xmax'])
                self.t.remove_row(i+1)
            else:
                i += 1


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

    def _cands_find(self, series, z_start, z_end, dz):
        trans = trans_parse(series)

        z_all = np.ravel([[to_z(x,t) for x in self.x] for t in trans])
        y_all = np.ravel([[self.y] for t in trans])

        (z_min, z_max) = (z_start, z_end) if z_start < z_end \
            else (z_end, z_start)
        z_int = z_all[np.logical_and(z_all>z_min, z_all<z_max)]
        y_int = y_all[np.logical_and(z_all>z_min, z_all<z_max)]

        if z_int == []: return z_int, y_int

        z_sort = np.sort(np.ravel(z_int))
        y_sort = y_int[np.argsort(np.ravel(z_int))]

        w_list = np.where(np.ediff1d(z_sort)<dz)[0]
        z_list = np.mean(np.vstack((z_sort[w_list], z_sort[w_list+1])), axis=0)
        y_list = np.mean(np.vstack((y_sort[w_list], y_sort[w_list+1])), axis=0)
        z_list = z_list if z_start<z_end else z_list[::-1]
        y_list = y_list if z_start<z_end else y_list[::-1]

        return z_list, y_list

    def _cands_find2(self, series, z_start, z_end, dz, single=False, logN=False):

        #print(z_start, z_end)
        # Compute all possible redshifts
        #trans = series_d[series]
        trans = trans_parse(series)
        if series == 'unknown':
            z_all = np.ravel([[self.x.to(au.nm)] for t in trans])
        else:
            z_all = np.ravel([[(x.to(au.nm)/xem_d[t].to(au.nm)).value-1. \
                                for x in self.x] for t in trans])
        trans_all = np.ravel([[t for x in self.x] for t in trans])
        y_all = np.ravel([[self.y] for t in trans])
        if logN:
            fosc_r = np.ravel([[fosc_d['Ly_a']/fosc_d[t]]*len(self.x) \
                               for t in trans])
            logN_all = np.ravel([[self.t['logN']] for t in trans]) \
                       + np.log10(fosc_r)
            #print('all')
            #print(z_all)
            #print(logN_all)

        # Select values within [z_start, z_end]
        (z_min, z_max) = (z_start, z_end) if z_start < z_end \
            else (z_end, z_start)
        z_sel = z_all[np.logical_and(z_all>z_min, z_all<z_max)]
        y_sel = y_all[np.logical_and(z_all>z_min, z_all<z_max)]
        trans_sel = trans_all[np.logical_and(z_all>z_min, z_all<z_max)]
        if logN:
            logN_sel = logN_all[np.logical_and(z_all>z_min, z_all<z_max)]
            #print('sel')
            #print(logN_sel)

        # Find coincidences
        #if len(series_d[series]) > 1:
        if len(trans) > 1:
            z_sort = np.sort(np.ravel(z_sel))
            y_sort = y_sel[np.argsort(np.ravel(z_sel))]
            trans_sort = trans_sel[np.argsort(np.ravel(z_sel))]
            if logN:
                logN_sort = logN_sel[np.argsort(np.ravel(z_sel))]
                #print('sort')
                #print(logN_sort)
            #print(z_sort)
            #print(trans_sort)

            w_range = np.where(np.logical_and(np.ediff1d(z_sort)<dz,
                                              trans_sort[:-1]!=trans_sort[1:]))[0]
            z_range = np.mean(np.vstack((z_sort[w_range], z_sort[w_range+1])),
                              axis=0)
            #y_range = y_sort[np.where(np.ediff1d(z_sort)<dz)]
            y_range = np.mean(np.vstack((y_sort[w_range], y_sort[w_range+1])),
                              axis=0)
            trans_range = trans_sort
            if logN:
                #logN_range = logN_sort[np.where(np.ediff1d(z_sort)<dz)]
                logN_range = np.log10(np.mean(np.vstack((
                    10**logN_sort[w_range], 10**logN_sort[w_range+1])),
                                              axis=0))
                #print('range 1')
                #print(z_range)
                #print(logN_range)

            z_range = z_range if z_start<z_end else z_range[::-1]
            y_range = y_range if z_start<z_end else y_range[::-1]
            trans_range = trans_range if z_start<z_end else trans_range[::-1]
            if logN:
                logN_range = logN_range if z_start<z_end else logN_range[::-1]
                #print('range 2')
                #print(z_range)
                #print(logN_range)
        else:
            z_range = z_sel
            y_range = y_sel
            trans_range = trans_sel
            if logN:
                logN_range = logN_sel

        if len(z_range) > 0:
            z_single = z_range[np.argmin(y_range)]
            trans_single = trans_range[np.argmin(y_range)]
            if logN:
                logN_single = logN_range[np.argmin(y_range)]
        else:
            z_single = None
            trans_single = None
            if logN:
                logN_single = None

        if single:
            if logN:
                return z_single, logN_single, trans_single
            else:
                return z_single, None, trans_single
        if logN:
            return z_range, logN_range, trans_range
        else:
            return z_range, [None]*len(z_range), trans_range


    def _syst_cand(self, series, z_start, z_end, dz, single=False, logN=False):

        # Compute all possible redshifts
        #trans = series_d[series]
        trans = trans_parse(series)
        if series == 'unknown':
            z_all = np.ravel([[self.x.to(au.nm)] for t in trans])
        else:
            z_all = np.ravel([[(x.to(au.nm)/xem_d[t].to(au.nm)).value-1. \
                                for x in self.x] for t in trans])
        y_all = np.ravel([[self.y] for t in trans])
        if logN:
            fosc_r = np.ravel([[fosc_d['Ly_a']/fosc_d[t]]*len(self.x) \
                               for t in trans])
            logN_all = np.ravel([[self.t['logN']] for t in trans]) \
                       + np.log10(fosc_r)
            #print('all')
            #print(z_all)
            #print(logN_all)

        # Select values within [z_start, z_end]
        (z_min, z_max) = (z_start, z_end) if z_start < z_end \
            else (z_end, z_start)
        z_sel = z_all[np.logical_and(z_all>z_min, z_all<z_max)]
        y_sel = y_all[np.logical_and(z_all>z_min, z_all<z_max)]
        if logN:
            logN_sel = logN_all[np.logical_and(z_all>z_min, z_all<z_max)]
            #print('sel')
            #print(logN_sel)

        # Find coincidences
        if len(series_d[series]) > 1:
            z_sort = np.sort(np.ravel(z_sel))
            y_sort = y_sel[np.argsort(np.ravel(z_sel))]
            if logN:
                logN_sort = logN_sel[np.argsort(np.ravel(z_sel))]
                #print('sort')
                #print(z_sort)
                #print(logN_sort)

            w_range = np.where(np.ediff1d(z_sort)<dz)[0]
            z_range = np.mean(np.vstack((z_sort[w_range], z_sort[w_range+1])),
                              axis=0)
            #y_range = y_sort[np.where(np.ediff1d(z_sort)<dz)]
            y_range = np.mean(np.vstack((y_sort[w_range], y_sort[w_range+1])),
                              axis=0)
            if logN:
                #logN_range = logN_sort[np.where(np.ediff1d(z_sort)<dz)]
                logN_range = np.log10(np.mean(np.vstack((
                    10**logN_sort[w_range], 10**logN_sort[w_range+1])),
                                              axis=0))
                #print('range 1')
                #print(z_range)
                #print(logN_range)

            z_range = z_range if z_start<z_end else z_range[::-1]
            y_range = y_range if z_start<z_end else y_range[::-1]
            if logN:
                logN_range = logN_range if z_start<z_end else logN_range[::-1]
                #print('range 2')
                #print(z_range)
                #print(logN_range)
        else:
            z_range = z_sel
            y_range = y_sel
            if logN:
                logN_range = logN_sel

        if len(z_range) > 0:
            z_single = z_range[np.argmin(y_range)]
            if logN:
                logN_single = logN_range[np.argmin(y_range)]
        else:
            z_single = None
            if logN:
                logN_single = None

        if single:
            if logN:
                return z_single, logN_single
            else:
                return z_single
        if logN:
            return z_range, logN_range
        else:
            return

    def append_replace(self, append, sess):
        if append and sess.lines is not None and len(sess.lines.t)>0:
            sess.lines._append(self)
            sess.lines._clean()
        else:
            sess.lines = self
