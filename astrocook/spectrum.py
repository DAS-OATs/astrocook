from .frame import Frame
from .functions import x_convert
from .line_list import LineList
#from .syst_list import SystList
from .message import *
from .vars import *
from astropy import units as au
from astropy.modeling.models import BlackBody
from astropy.modeling.powerlaws import PowerLaw1D
from astropy.table import Column
from astropy.stats import sigma_clip
import bisect
#from astropy import constants as aconst
#from astropy import table as at
from copy import deepcopy as dc
from time import time
import logging
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import argrelmin, argrelmax, fftconvolve
from scipy.interpolate import UnivariateSpline as uspline
from scipy.stats import sem
from tqdm import tqdm

class Spectrum(Frame):
    """Class for spectra

    A Spectrum is a Frame with methods for handling spectral operations."""

    def __init__(self,
                 x=[],
                 xmin=[],
                 xmax=[],
                 y=[],
                 dy=[],
                 xunit=au.nm,
                 yunit=au.erg/au.cm**2/au.s/au.nm,
                 meta={},
                 dtype=float,
                 cont=[],
                 resol=[]):
        super(Spectrum, self).__init__(x, xmin, xmax, y, dy, xunit, yunit, meta,
                                       dtype)
        if cont != []:
            self._t['cont'] = cont*self._yunit
        if resol != []:
            self._t['resol'] = resol


    def _copy(self, sel=None):
        copy = super(Spectrum, self)._copy(sel)
        cols = [c for c in self._t.colnames \
                if c not in ['x', 'xmin', 'xmax', 'y', 'dy']]
        for c in cols:
            copy._t[c] = self._t[c][sel]
        return copy

    def _deredden(self, ebv=0.03, rv=3.1):

        invx = 1/self.x.to(au.micron).value
        a = np.zeros(len(invx))
        b = np.zeros(len(invx))

        # IR
        ir_w = np.where(np.logical_and(invx>0.3, invx<1.1))
        a[ir_w] = 0.574 * invx[ir_w]**1.61
        b[ir_w] = -0.527 * invx[ir_w]**1.61

        # Visual/NIR (0'Donnell 1994)
        vis_w = np.where(np.logical_and(invx>1.1, invx<3.3))
        c1 = [1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
        c2 = [0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]
        a[vis_w] = np.polyval(c1[::-1], invx[vis_w]-1.82)
        b[vis_w] = np.polyval(c2[::-1], invx[vis_w]-1.82)

        # Mid UV
        muv_w = np.where(np.logical_and(invx>3.3, invx<8.0))
        f_a = np.zeros(len(muv_w[0]))
        f_b = np.zeros(len(muv_w[0]))
        f_w = np.where(invx[muv_w]>5.9)
        invx_w = invx[muv_w][f_w]-5.9
        f_a[f_w] = -0.04473 * invx_w**2 - 0.009779 * invx_w**3
        f_b[f_w] = 0.2130 * invx_w**2 + 0.1207 * invx_w**3
        a[muv_w] = 1.752 - 0.316*invx[muv_w] \
                   - (0.104 / ((invx[muv_w]-4.67)**2 + 0.341)) + f_a
        b[muv_w] = -3.090 + 1.825*invx[muv_w] \
                   + (1.206 / ((invx[muv_w]-4.62)**2 + 0.263)) + f_b

        # Far UV
        fuv_w = np.where(np.logical_and(invx>8.0, invx<11.0))
        c1 = [-1.073, -0.628, 0.137, -0.070]
        c2 = [13.670, 4.257, -0.420, 0.374]
        a[fuv_w] = np.polyval(c1[::-1], invx[fuv_w]-8.0)
        b[fuv_w] = np.polyval(c2[::-1], invx[fuv_w]-8.0)

        av = rv*ebv
        al = av * (a+b/rv)

        self.y = self.y * 10**(0.4*al)

        return 0


    def _flux_ccf(self, col1, col2, dcol1, dcol2, vstart, vend, dv,
                  weighted=False):
        vstart = vstart.to(au.km/au.s).value
        vend = vend.to(au.km/au.s).value
        dv = dv.to(au.km/au.s).value
        sd = -1*int(np.floor(np.log10(dv)))-1
        spec_x = self.x.value[:]
        #spec_x_2 = self.x.value[:]*(1+0.063/aconst.c.to(au.km/au.s).value)

        xmin = spec_x[~np.isnan(spec_x)][0]
        xmax = spec_x[~np.isnan(spec_x)][-1]
        dv_orig = (self._t['xmax']-self._t['xmin'])/spec_x*aconst.c.to(au.km/au.s).value
        xmean = 0.5*(xmin+xmax)
        v_shift = np.arange(vstart, vend+dv, dv)

        x_shift = xmean * v_shift/aconst.c.to(au.km/au.s).value
        xstart = xmean * vstart/aconst.c.to(au.km/au.s).value
        xend = xmean * vend/aconst.c.to(au.km/au.s).value
        dx = xmean * dv/aconst.c.to(au.km/au.s).value

        scale = int(np.rint(np.nanmedian(dv_orig)/dv))

        x_osampl = np.arange(xmin+xstart, xmax+xend, dx)
        y1_osampl = np.interp(x_osampl, spec_x, self._t[col1])
        y2_osampl = np.interp(x_osampl, spec_x, self._t[col2])
        #y2_osampl = np.interp(x_osampl, spec_x_2, self._t[col2])
        dy1_osampl = np.interp(x_osampl, spec_x, self._t[dcol1])
        dy2_osampl = np.interp(x_osampl, spec_x, self._t[dcol2])
        #dy2_osampl = np.interp(x_osampl, spec_x_2, self._t[dcol2])

        pan = len(x_shift)//2
        pan_l, pan_r = int(abs(len(x_shift)*xstart/np.abs(xend-xstart))), \
            int(abs(len(x_shift)*xend/np.abs(xend-xstart)))

        #print(scale, pan, pan_l, pan_r)
        ccf = []
        chi2 = []
        chi2r = []
        check = np.inf
        for i, xs in enumerate(x_shift):
            """
            y1 = y1_osampl[pan:-pan-1]-np.nanmean(y1_osampl)
            y2 = y2_osampl[i:-2*pan+i-1]-np.nanmean(y2_osampl)
            dy = dy1_osampl[pan:-pan-1]

            ccf.append(np.nanmean(y2 * y1)/np.sqrt(np.nanmean(y2**2) * np.nanmean(y1**2)))
            """

            x = x_osampl+xs

            y1 = y1_osampl[pan_l:-pan_r-1]
            y2 = y2_osampl[i:-pan_l-pan_r+i-1]

            dy = dy1_osampl[pan_l:-pan_r-1]

            y1 = y1[::scale]
            y2 = y2[::scale]
            dy = dy[::scale]

            y1m = y1-np.nanmedian(y1)
            y2m = y2-np.nanmedian(y2)

            ccf.append(np.nanmean(y2m*y1m)/np.sqrt(np.nanmean(y2m**2) * np.nanmean(y1m**2)))
            #"""
            chi2i = (y1-y2)**2/dy**2

            if weighted:
                bf = np.abs(np.gradient(y2))
                bf = bf * len(chi2i)/np.sum(bf)
                chi2i = chi2i * bf

            chi2i_sum = np.nansum(chi2i)

            chi2_plot = False
            if chi2_plot:
                chi2ran = np.arange(0,50,0.5)
                from scipy.stats import chi2 as scipychi2
                #max_chi2 = np.argmin(np.abs(scipychi2.pdf(chi2ran, 1)-1/len(chi2i)))
                plt.hist(chi2i, bins=chi2ran, density=True)
                plt.plot(chi2ran, scipychi2.pdf(chi2ran, 1), color='red')
                plt.xlim(0,20)
                plt.ylim(1e-5,1)
                plt.yscale('log')
                plt.xlabel(r'$\chi_r^2$')
                plt.ylabel('Frequency')
                plt.show()

            resid_plot = False
            if resid_plot:
                if chi2i_sum < check:
                    check = chi2i_sum
                else:
                    sss = chi2i>15
                    """
                    from astropy.table import Table
                    t = Table([x[pan_l:-pan_r-1][::scale][sss]])
                    t['g'] = np.floor(100*t['col0'])/100
                    t['g2'] = np.floor(100*t['col0'])/100+0.01
                    #t.pprint(max_lines=375)
                    t2 = Table([np.unique(t['g'])])
                    t2['t'] = np.floor(100*t2['g']+1)/100
                    rem = []
                    for i in range(len(t2)-1,0,-1):
                        if t2['g'][i]-t2['g'][i-1]<0.015:
                            t2['t'][i-1]=t2['t'][i]
                            rem.append(i)
                    t2.remove_rows(rem)
                    t2.pprint(max_lines=1000)
                    """
                    plt.scatter(x[pan_l:-pan_r-1][::scale][sss], chi2i[sss], s=1, color='black')
                    plt.show()

            chi2.append(chi2i_sum)
            chi2r.append(chi2i_sum/len(y1))

        return np.array(v_shift), np.array(ccf), np.array(chi2), np.array(chi2r)


    def _gauss_convolve(self, std=20, input_col='y', output_col='conv',
                        verb=True):

        # Create profile
        self._x_convert()
        x = self._safe(self.x)
        mean = np.median(x)
        prof = np.exp(-((x - mean) / std).value**2)
        if (len(prof) % 2 == 0):
            prof = prof[:-1]
        prof = prof / np.sum(prof)

        # Convolve
        if verb:
            if output_col not in self._t.colnames:
                logging.info("I'm adding column '%s'." % output_col)
        conv = dc(self._t[input_col])

        # Add padding to avoid edge issues
        len_conv = len(conv)
        conv = Column(np.append(np.ones(len_conv)*conv[0], conv))
        conv = Column(np.append(conv, np.ones(len_conv)*conv[-1]))

        safe = np.array(self._safe(conv), dtype=float)
        mode = 'same'
        try:
            conv[self._where_safe] = fftconvolve(safe, prof, mode=mode)\
                                                 *self._t[input_col].unit
        except:
            conv[self._where_safe] = fftconvolve(safe, prof, mode=mode)
        self._t[output_col] = conv[len_conv:-len_conv] * self._t[input_col].unit
        self._x_convert(xunit=self._xunit_old)
        self._xunit_old = self._xunit

        return 0


    """
    def _syst_fit(self, series='CIV', z=1.6971, logN=13, b=10, resol=45000):
        @brief Add a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @return System list

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)

        self._systs = SystList(self)
        self._systs._add_fit(series, z, logN, b, resol)

        return self._systs
    """


    def _lines_mask(self, lines, source=None):
        """ @brief Create a mask consisting on the ['xmin', 'xmax'] regions from
        the associated line list
        @return 0
        """

        x = self._safe(self.x)
        mask = np.zeros(len(x), dtype=bool)

        if source is not None:
            where = lines.t['source'] == source
        else:
            where = range(len(lines.t))

        lines_xmin = lines.xmin[where]
        lines_xmax = lines.xmax[where]

        for (xmin, xmax) in zip(lines_xmin, lines_xmax):
            mask += np.logical_and(x>=xmin, x<=xmax)
        if 'lines_mask' not in self._t.colnames:
            logging.info("I'm adding column 'lines_mask' to spectrum.")
            self._t['lines_mask'] = np.empty(len(self.x), dtype=bool)
        self._t['lines_mask'][self._where_safe] = mask

        return 0


    def _node_add(self, nodes, x, y):
        if isinstance(nodes, str):
            spl = nodes.split('.')
            nodes = getattr(getattr(self._gui, spl[0]), spl[1])
        abs = np.abs(self.x.to(self._xunit).value-x)
        sel = abs[~np.isnan(abs)].argmin()
        row = []
        for c in nodes.t.colnames:
            row.append(y) if c == 'y' else row.append(self.t[sel][c])
        nodes.t.add_row(row)
        nodes.t.sort('x')


    def _node_remove(self, nodes, x):
        if isinstance(nodes, str):
            spl = nodes.split('.')
            nodes = getattr(getattr(self._gui, spl[0]), spl[1])
        sel = np.abs(nodes.x.to(self._xunit).value-x).argmin()
        nodes.t.remove_rows(sel)


    def _nodes_clean(self, nodes, kappa=5.0):
        """
        y_sel = [True]
        hw = window//2
        while np.any(y_sel):
            y_start = nodes.y.value
            yg = [y_start[max(0, i-hw):i+hw+1] for i in range(len(y_start))]
            yg_median = [np.median(g) for g in yg]
            yg_std = [np.std(g) for g in yg]

            y_sel = [np.abs(y-m) > kappa*s \
                     for y, m, s in zip(y_start, yg_median, yg_std)]
            #print(yg)
            #print(yg_median)
            #print(yg_std)
            #print(y_sel)
            print(np.sum(y_sel))
            #print(nodes._t[y_sel])
            nodes._t.remove_rows(y_sel)
        """
        lines = nodes._peaks_find(col='y', kappa=kappa, kind='max')
        y_sel = np.where(np.in1d(nodes.x.value, lines.x.value))[0]
        nodes._t.remove_rows(y_sel)
        while np.any(y_sel):
            lines = nodes._peaks_find(col='y', kappa=kappa)
            y_sel = np.where(np.in1d(nodes.x.value, lines.x.value))[0]
            #print(len(y_sel))
            nodes._t.remove_rows(y_sel)
        #"""
        return nodes


    def _nodes_extract(self, delta_x=1500, xunit=au.km/au.s, mode='std'):

        self._slice(delta_x, xunit)
        x_ave = []
        xmin_ave = []
        xmax_ave = []
        y_ave = []
        dy_ave = []
        if 'lines_mask' not in self._t.colnames:
            logging.warning("Lines weren't masked. I'm taking all spectrum.")


        for s in self._slice_range:
            if mode=='std':
                try:
                    where_s = np.where(np.logical_and(self._t['slice']==s,
                                                    self._t['lines_mask']==0))
                except:
                    where_s = np.where(self._t['slice']==s)
            elif mode=='cont':
                where_s = np.where(self._t['slice']==s)

            # Use deabs column if present
            if mode == 'std':
                y = self._t['deabs'] if 'deabs' in self._t.colnames else self.y.value
            elif mode == 'cont':
                y = self._t['cont']


            if len(where_s[0])>0.1*len(np.where(self._t['slice']==s)[0]):
                x_where_s = self.x[where_s].value
                y_where_s = y[where_s]
                dy_where_s = self.dy[where_s].value
                x_ave.append(np.median(x_where_s))
                xmin_ave.append(x_where_s[0])
                xmax_ave.append(x_where_s[-1])
                if mode == 'std':
                    y_ave.append(np.median(y_where_s))
                elif mode == 'cont':
                    y_ave.append(np.interp(np.median(x_where_s), x_where_s, y_where_s))
                dy_ave.append(sem(y_where_s))
        x = np.array(x_ave) * self._xunit
        xmin = np.array(xmin_ave) * self._xunit
        xmax = np.array(xmax_ave) * self._xunit
        y = np.array(y_ave) * self._yunit
        dy = np.array(dy_ave) * self._yunit


        return Spectrum(x, xmin, xmax, y, dy, self._xunit, self._yunit)


    def _nodes_interp(self, lines, nodes, smooth=0):
        """ @brief Interpolate nodes with a univariate spline to estimate the
        emission level.
        @param smooth Smoothing of the spline
        @return 0
        """

        if isinstance(lines, str):
            spl = nodes.split('.')
            lines = getattr(getattr(self._gui, spl[0]), spl[1])
        if isinstance(nodes, str):
            spl = nodes.split('.')
            nodes = getattr(getattr(self._gui, spl[0]), spl[1])

        x = nodes.x.value
        y = nodes.y.value
        dy = nodes.dy.value
        #isnan = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),
        #                      np.isnan(dy))
        isnan = np.logical_or(np.isnan(x),np.isnan(y))
        dy[np.isnan(dy)] = np.median(dy[~np.isnan(dy)])
        if np.sum(np.isnan(dy)) > 0:
            spl = uspline(x[~isnan], y[~isnan], s=smooth)
        else:
            spl = uspline(x[~isnan], y[~isnan], w=dy[~isnan], s=smooth)
        cont = spl(self.x)*self._yunit
        logging.info("I'm using interpolation as continuum.")
        if 'cont' not in self._t.colnames:
            logging.info("I'm adding column 'cont'.")
        self._t['cont'] = cont #spl(self.x)
        if hasattr(lines, '_t'):
            lines._t['cont'] = np.interp(lines.x, self.x, cont)
        nodes._t['cont'] = np.interp(nodes.x, self.x, cont)
        return 0


    def _peaks_find(self, col='conv', kind='min', kappa=3.0, **kwargs):
        y = self._safe(self._t[col])
        min_idx = np.hstack(argrelmin(y, **kwargs))
        max_idx = np.hstack(argrelmax(y, **kwargs))
        ext_idx = np.sort(np.append(min_idx, max_idx))
        ext = self._copy(ext_idx)
        if len(ext.t) > 0:
            # Set xmin and xmax from adjacent extrema
            ext.xmin = np.append(self.x[0], ext.x[:-1])
            ext.xmax = np.append(ext.x[1:], self.x[-1])

            diff_y_left = (ext._t[col][:-2] - ext._t[col][1:-1])
            diff_y_right = (ext._t[col][2:] - ext._t[col][1:-1])
            if kind == 'max':
                diff_y_left = -diff_y_left
                diff_y_right = -diff_y_right

        # Check if the difference is above threshold
        #for m,M,l,r in zip(ext.xmin, ext.xmax, diff_y_left, diff_y_right):

        #    print(m,M,l,r)
            diff_y_max = np.minimum(diff_y_left, diff_y_right)

        # +1 is needed because sel is referred to the [1:-1] range of rows
        # in the spectrum
            dy = ext.dy[1:-1]
            dy[np.isnan(dy)] = np.median(dy[~np.isnan(dy)])
            sel = np.where(np.greater(diff_y_max, dy * kappa))[0]+1
        else:
            sel = []

        lines = ext._copy(sel)

        return lines

    def _rebin(self, xstart, xend, dx, xunit, y, dy, kappa=5, filling=np.nan):

        # Convert spectrum into chosen unit
        # A deep copy is created, so the original spectrum is preserved


        self.t.sort('x')
        self._x_convert(xunit=xunit)

        # Create x, xmin, and xmax
        from .format import Format
        format = Format()
        if xstart is None:
            xstart = np.nanmin(self.x)
        else:
            xstart = (xstart*au.nm).to(xunit, equivalencies=equiv_w_v)
        if xend is None:
            xend = np.nanmax(self.x)
        else:
            xend = (xend*au.nm).to(xunit, equivalencies=equiv_w_v)
        x = np.arange(xstart.value, xend.value, dx) * xunit
        xmin, xmax = format._create_xmin_xmax(x)

        # Compute y and dy combining contributions
        im = 0
        iM = 1
        xmin_in = self.xmin[iM].value
        xmax_in = self.xmax[im].value
        y_out = np.array([]) * y.unit
        dy_out = np.array([]) * y.unit
        print_time = False
        xmin_value = np.array(self.xmin.value)
        xmax_value = np.array(self.xmax.value)
        for i, (m, M) \
            in enum_tqdm(zip(xmin.value, xmax.value), len(xmin),
                         "spectrum: Rebinning"):
            if print_time:
                print('')
                t1 = time()
                print(t1)
            """
            while xmin_in < M:
                iM += 1
                try:
                    xmin_in = self.xmin[iM].value
                    #print(xmin_in)
                except:
                    break
            while xmax_in < m:
                im += 1
                try:
                    xmax_in = self.xmax[im].value
                    #print(xmax_in)
                except:
                    break
            """
            im = bisect.bisect_left(xmax_value, m)
            iM = bisect.bisect_right(xmin_value, M)
            #print('im  ',im, iM)
            if print_time:
                t15 = time()
                print(t15, t15-t1)

            frac = (np.minimum(M, xmax_value[im:iM])\
                    -np.maximum(m, xmin_value[im:iM]))/dx
            if print_time:
                t17 = time()
                print(t17, t17-t15)
            ysel = y[im:iM]
            #print(m, M, self.xmin[im:iM], self.xmax[im:iM])
            #print(frac)

            #print(frac[w],frac)
            dysel = dy[im:iM]

            nw = np.where(~np.isnan(ysel))
            ysel = ysel[nw]
            dysel = dysel[nw]
            frac = frac[nw]
            #print(dysel)
            #mask = sigma_clip(ysel, masked=True).mask
            #if np.sum(~mask)>0:
            #    frac = frac[~mask]
            #    ysel = ysel[~mask]
            #    dysel = dysel[~mask]
            from astropy.stats import sigma_clip

            # Optional kappa-sigma clipping of outliers
            if kappa is not None:
                yclip = sigma_clip(ysel, sigma=kappa, masked=True)
                if len(yclip)>0:
                    w = np.where(np.logical_and(frac>0, yclip.mask==False))
                else:
                    w = np.where(frac>0)
            else:
                w = np.where(frac>0)

            if print_time:
                t2 = time()
                print(t2, t2-t16)

            if len(frac[w]) > 0:
                weights = (frac[w]/dysel[w]**2).value
                #print(frac[w], np.sum(frac[w])/len(frac[w]))
                #nw = np.where(~np.isnan(ysel[w]))
                #print(weights)
                #print(frac[w])
                if np.any(np.isnan(dysel)) or np.any(dysel==0.0) or np.sum(weights)==0.0:# and False:
                    y_out = np.append(y_out, np.average(ysel[w], weights=frac[w]))
                else:
                    y_out = np.append(y_out, np.average(ysel[w], weights=weights))
                    #y_out = np.append(y_out, np.average(ysel[w], weights=frac[w]/dysel[w]**2))
                dy_out = np.append(dy_out, np.sqrt(np.nansum(weights**2*dysel[w].value**2))\
                                                   /np.nansum(weights)*y.unit)
                #dy_out = np.append(dy_out, np.sqrt(np.sum(frac**2/dysel**2))\
                #                                   /np.sum(frac/dysel**2))
            else:
                y_out = np.append(y_out, filling)
                dy_out = np.append(dy_out, filling)
            if print_time:
                t3 = time()
                print(t3, t3-t2)

        # Create a new spectrum and convert it to the units of the original one
        out = Spectrum(x, xmin, xmax, y_out, dy_out, xunit=xunit, yunit=y.unit,
                       meta=self.meta)
        out._x_convert(xunit=self._xunit_old)
        self._x_convert(xunit=self._xunit_old)
        self._xunit_old = self._xunit
        return out


    def _resol_est(self, px, update):
        dx = self.xmax[px-1:].value-self.xmin[:1-px].value
        if px%2:
            dx = np.append(np.append([dx[0]]*(px//2), dx), [dx[-1]*(px//2)])\
                 *self.x.unit
        else:
            dx = np.append([dx[0]]*(px//2), dx)*self.x.unit

        if update:
            if 'resol' not in self._t.colnames:
                logging.info("I'm adding column 'resol'.")
            else:
                logging.info("I'm updating column 'resol'.")
            self._t['resol'] = self.x/dx
        logging.info('The mean estimated resolution is %4.2f.'
                     % np.nanmean(self._t['resol']))
        return 0


    def _slice(self, delta_x=1000, xunit=au.km/au.s):
        """ @brief Create 'slice' columns. 'slice' columns contains an
        increasing counter to split 'x' values into evenly-sized slices
        (typically defined in velocity space).
        @param delta_x Size of slices
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        self._x_convert(xunit=xunit)
        x = self._safe(self.x)
        self._t['slice'] = np.empty(len(self.x), dtype=int)
        self._t['slice'][self._where_safe] = np.array(x//delta_x)
        self._slice_range = range(self._t['slice'][self._where_safe][0],
                                  self._t['slice'][self._where_safe][-1])
        self._x_convert(xunit=self._xunit_old)
        self._xunit_old = self._xunit
        return 0


    def _stats_print(self, xmin=0, xmax=np.inf):
        try:
            xmin = xmin.to(au.nm).value
            xmax = xmax.to(au.nm).value
        except:
            pass
        sel = np.where(np.logical_and(self.x.to(au.nm).value > xmin,
                                      self.x.to(au.nm).value < xmax))
        x = self.x[sel]
        dx = self.xmax[sel]-self.xmin[sel]
        dxv = x_convert(self.xmax[sel], xunit=au.km/au.s) \
              -x_convert(self.xmin[sel], xunit=au.km/au.s)
        y = self.y[sel]
        dy = self.dy[sel]
        self._stats = {'min_x': np.nanmin(x),
                       'max_x': np.nanmax(x),
                       'mean_x': np.nanmean(x),
                       'min_dx': np.nanmin(dx),
                       'max_dx': np.nanmax(dx),
                       'mean_dx': np.nanmean(dx),
                       'min_dxv': np.nanmin(dxv),
                       'max_dxv': np.nanmax(dxv),
                       'mean_dxv': np.nanmean(dxv),
                       'min_y': np.nanmin(y),
                       'max_y': np.nanmax(y),
                       'mean_y': np.nanmean(y),
                       'median_y': np.nanmedian(y),
                       'std_y': np.nanstd(y),
                       'min_dy': np.nanmin(dy),
                       'max_dy': np.nanmax(dy),
                       'mean_dy': np.nanmean(dy),
                       'median_dy': np.nanmedian(dy.value)*dy.unit}
        #self._stats_tup = tuple(np.ravel([(self._stats[s].value,
        #                                  self._stats[s].unit) \
        #                                  for s in self._stats]))
        self._stats_tup = (self._stats['min_x'].value,
                           self._stats['max_x'].value,
                           self._stats['min_x'].unit,
                           self._stats['mean_x'].value,
                           self._stats['mean_x'].unit,
                           self._stats['min_dx'].value,
                           self._stats['max_dx'].value,
                           self._stats['min_dx'].unit,
                           self._stats['min_dxv'].value,
                           self._stats['max_dxv'].value,
                           self._stats['min_dxv'].unit,
                           self._stats['mean_dx'].value,
                           self._stats['mean_dx'].unit,
                           self._stats['mean_dxv'].value,
                           self._stats['mean_dxv'].unit,
                           self._stats['min_y'].value,
                           self._stats['max_y'].value,
                           self._stats['min_y'].unit,
                           self._stats['mean_y'].value,
                           self._stats['mean_y'].unit,
                           self._stats['median_y'].value,
                           self._stats['median_y'].unit,
                           self._stats['std_y'].value,
                           self._stats['std_y'].unit,
                           self._stats['min_dy'].value,
                           self._stats['max_dy'].value,
                           self._stats['min_dy'].unit,
                           self._stats['mean_dy'].value,
                           self._stats['mean_dy'].unit)
        self._stats_text = "Statistics in the selected region:\n" \
                           " x range:        %3.4f - %3.4f %s\n" \
                           " x mean:         %3.4f %s\n" \
                           " bin size range: %3.4f - %3.4f %s (%3.4f - %3.4f %s)\n" \
                           " bin size mean:  %3.4f %s (%3.4f %s)\n" \
                           " y range:        %3.4e - %3.4e %s\n" \
                           " y mean:         %3.4e %s\n" \
                           " y median:       %3.4e %s\n" \
                           " y stdev:        %3.4e %s\n" \
                           " dy range:       %3.4e - %3.4e %s\n" \
                           " dy mean:        %3.4e %s" \
                           % self._stats_tup
        if xmin == 0.0 and np.isinf(xmax):
            region = ("ALL SPECTRUM",)
        else:
            temp = (self._stats['min_x'].value, self._stats['max_x'].value,
                    self._stats['min_x'].unit)
            region = ("REGION: %3.4f-%3.4f %s" % temp,)
        red = region+(self._stats['mean_y'].value,
                      self._stats['mean_y'].unit,
                      self._stats['median_y'].value,
                      self._stats['median_y'].unit,
                      self._stats['std_y'].value,
                      self._stats['std_y'].unit)
        self._stats_text_red = "%s\n" \
                               "y mean: %3.4e %s\n" \
                               "y median: %3.4e %s\n" \
                               "y stdev: %3.4e %s" \
                               % red
        for s in self._stats_text.split('\n'):
            logging.info(s)


    def _template_bb(self, temp=6000, scale=1.0):
        bb = BlackBody(temperature=temp*au.K, scale=scale*au.erg/(self._xunit*au.cm**2*au.s*au.sr))
        output_col = 'blackbody'
        if output_col not in self._t.colnames:
            logging.info("I'm adding column '%s'." % output_col)
        self._t[output_col] = bb(self.x)


    def _template_pl(self, ampl=1.0, x_ref=None, index=-1.0):
        if x_ref == None: x_ref = np.mean(self.x)
        pl = PowerLaw1D(amplitude=ampl, x_0=x_ref, alpha=-index)
        output_col = 'power_law'
        if output_col not in self._t.colnames:
            logging.info("I'm adding column '%s'." % output_col)
        self._t[output_col] = pl(self.x)


    def _y_scale(self, fact):
        super(Spectrum, self)._y_scale(fact)

        cols = ['model', 'deabs', 'y_rm', 'cont']
        for c in cols:
            if c in self._t.colnames:
                self._t[c] = self._t[c] * fact
        return 0


    def _zap(self, xmin, xmax):

        xmin = np.ravel(np.array(xmin))
        if xmax is not None:
            xmax = np.ravel(np.array(xmax))
            for m, M in zip(xmin, xmax):
                w = np.where(np.logical_and(self.x.value>m, self.x.value<M))
                self._t['y'][w] = np.interp(
                                      self.x[w].value,
                                      [self.x[w][0].value, self.x[w][-1].value],
                                      [self.y[w][0].value, self.y[w][-1].value])*self._yunit
                self._t['dy'][w] = np.interp(
                                      self.x[w].value,
                                      [self.x[w][0].value, self.x[w][-1].value],
                                      [self.dy[w][0].value, self.dy[w][-1].value])*self._yunit
        else:
            #self._t.remove_row(np.abs(self._t['x'] - xmin).argmin())
            for x in xmin:
                r = np.nanargmin(np.abs(self._t['x'] - x))
                self._t['x'][r] = np.nan

        return 0
