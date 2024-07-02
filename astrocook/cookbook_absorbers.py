from .vars import resol_def
from .line_list import LineList

import logging
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import argrelmin, argrelmax, find_peaks
from scipy.special import erf, erfc, erfinv
import sys
import time

prefix = "[INFO] cookbook_absorbers:"

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


class CookbookAbsorbers(object):
    """ Cookbook of utilities for modeling absorbers
    """

    def __init__(self):
        super(CookbookAbsorbers, self).__init__()
        self._refit_n = 3
        self._chi2rav_thres = 1e-2
        self._chi2r_thres = np.inf
        self._dlogN_thres = np.inf
        self._max_nfev = max_nfev_def
        self._sel_fit = False

    def _lines_cands_find(self, series, z_start, z_end, dz):
        return self.sess.lines._cands_find(series, z_start, z_end, dz)


    def lines(self, kind='abs', prominence=None, append=True):
        """ @brief Find lines

        @details Find absorption or emission lines, based on their prominence.

        @param kind Kind of lines (`abs` or `em`)
        @param prominence Prominence of lines (as in `scipy.signal.find_peaks`)
        @param append Append lines to existing line list
        @return 0
        """

        try:
            kind = str(kind)
            prominence = None if prominence in [None, 'None'] else float(prominence)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if kind not in ['abs', 'em']:
            logging.error("`kind` should be `abs` or `em`. Aborting.")
            return 0

        spec = self.sess.spec
        fact = -1 if kind=='abs' else 1

        ynorm = fact*(spec._t['y'])
        if prominence is None:
            prominence = 5*(spec._t['dy'])
        peaks, properties = find_peaks(ynorm, prominence=prominence)
        lines = LineList(spec._t['x'][peaks],
                         spec._t['xmin'][peaks],
                         spec._t['xmax'][peaks],
                         spec._t['y'][peaks],
                         spec._t['dy'][peaks],
                         spec._t['cont'][peaks],
                         ['y']*len(peaks),
                         spec._xunit, spec._yunit,
                         meta=spec._meta)
        if append and self.sess.lines is not None \
            and len(self.sess.lines.t) > 0:
            self.sess.lines._append(lines)
            self.sess.lines._clean()
        else:
            logging.warning("I'm updating column 'like' in spectrum.")
        for c in spec._t.colnames:
            if 'like_' in c:
                spec._t['like'] = np.fmax(spec._t['like'], spec._t[c])


        return likes, z_likes

    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def systs_complete_from_like(self, series='all', series_ref=None, z_start=0,
                                 z_end=6, binz=1e-2, dz=1e-4,
                                 modul=10, sigma=2, distance=3,
                                 logN=logN_def, b=b_def, resol=resol_def,
                                 chi2r_thres=np.inf, dlogN_thres=np.inf,
                                 refit_n=0, chi2rav_thres=1e-2,
                                 max_nfev=max_nfev_def, append=True):
        """ @brief Complete systems from likelihood
        @details TBD
        @param series Series of transitions
        @param series_ref Reference series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param binz Bin size to group existing redshifts
        @param dz Threshold for redshift coincidence
        @param modul Modulation of the error function
        @param sigma Threshold for accepting systems
        @param distance Distance between systems in pixels
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @param append Append systems to existing system list
        @return 0
        """

        try:
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            binz = float(binz)
            dz = float(dz)
            modul = float(modul)
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = float(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0


        spec = self.sess.spec
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()

        if series == 'all':
            trans_ex = np.unique(np.ravel([trans_parse(s)
                                           for s in systs._t['series']]))
            trans_n = list(set(trans_d)-set(trans_ex))
            #print(trans_n)
            series = ';'.join(trans_n)

        if series_ref != None:
            w = np.where(systs._t['series']==series_ref)
            #hist, edges = np.histogram(systs._t['z'][w], bins=np.arange(0, 10, binz))

            # Rebecca's fix
            median_z = np.nanmedian(systs._t['z'][w])
            z_lower_boundary = np.round(((median_z + 0.5 * binz) % binz) / dz) * dz
            hist, edges = np.histogram(systs._t['z'][w], bins=np.arange(z_lower_boundary, 10, binz))
        else:
            hist, edges = np.histogram(systs._t['z'], bins=np.arange(0, 10, binz))

        w_z = np.logical_and(edges>z_start, edges<z_end)


        self._likes, self._z_likes = {}, {}
        for z in edges[:-1][hist>0]:
            if z>z_start and z+binz<z_end:
                z_s = z
                z_e = z+binz
            elif z<z_start and z+binz>z_end:
                z_s = z_start
                z_e = z_end
            elif z<z_start and z+binz>z_start:
                z_s = z_start
                z_e = z+binz
            elif z<z_end and z+binz>z_end:
                z_s = z
                z_e = z_end
            else:
                z_s = np.nan
                z_e = np.nan
            #print(z, z+binz, z_start, z_end, z_s, z_e)
            if not np.isnan(z_s) and not np.isnan(z_e):
                likes, z_likes = self._abs_like(series, 'y', z_s, z_e, dz, modul)
                """
            for s in likes.keys():
                if s not in self._likes.keys():
                    self._likes[s] = likes[s]
                    self._z_likes[s] = z_likes[s]
                else:
                    self._likes[s] = np.append(self._likes[s], likes[s])
                    self._z_likes[s] = np.append(self._z_likes[s], likes[s])
                """
                self._likes = likes
                self._z_likes = z_likes
                self._systs_like(series, sigma, distance, logN, b, resol, chi2r_thres,
                                 dlogN_thres, refit_n, chi2rav_thres, max_nfev, append)

        if compressed:
            systs._compress()
        """
        logging.info("I'm completing systems with %i additional transitions." \
                     % len(self._likes))
        self._systs_like(series, thres, logN, b, resol, chi2r_thres,
                         dlogN_thres, refit_n, chi2rav_thres, max_nfev, append)
        """
        return 0


    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def systs_new_from_like(self, series='Ly-a', col='y', z_start=0, z_end=6,
                            dz=1e-4, modul=10, sigma=2, distance=3,
                            logN=logN_def, b=b_def, resol=resol_def,
                            chi2r_thres=np.inf, dlogN_thres=np.inf,
                            refit_n=0, chi2rav_thres=1e-2, max_nfev=max_nfev_def,
                            append=True):

        """ @brief New systems from likelihood
        @details Add Voigt models by testing transitions in a given redshift
        range and assignign the most likely identification to absorption
        features.
        @param series Series of transitions
        @param col Column to apply the likelihood
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param modul Modulation of the error function
        @param sigma Threshold for accepting systems
        @param distance Distance between systems in pixels
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @param append Append systems to existing system list
        @return 0
        """

        try:
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            modul = float(modul)
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = float(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0

        self._likes, self._z_likes = self._abs_like(series, col, z_start, z_end, dz,
                                                    modul)
        self._systs_like(series, sigma, distance, logN, b, resol, chi2r_thres,
                         dlogN_thres, refit_n, chi2rav_thres, max_nfev, append)

        return 0


    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def _systs_like(self, series='Ly-a', sigma=2, distance=3, logN=logN_def,
                    b=b_def, resol=resol_def, chi2r_thres=np.inf,
                    dlogN_thres=np.inf, refit_n=0, chi2rav_thres=1e-2,
                    max_nfev=max_nfev_def, append=True):

        try:
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = float(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if series == 'all':
            series = ';'.join(trans_d)

        spec = self.sess.spec
        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()


        likes = self._likes
        z_likes = self._z_likes
        series_split = series.split(';')

        def _nan_check(spec, s_list, z_list, logN_list, resol_list):
            # Check that components do not fall in masked regions
            sel = []
            for ssel,zsel in zip(s_list, z_list):
                ysel = []
                for t in trans_parse(ssel):
                    xsel = to_x(zsel, t)
                    ysel.append(np.interp(xsel, spec._t['x'],
                                          spec._t['y']))
                sel.append(not np.any(np.isnan(ysel)))
            wsel = np.where(sel)[0]
            s_list = np.array(s_list)[wsel]
            z_list = np.array(z_list)[wsel]
            logN_list = np.array(logN_list)[wsel]
            resol_list = np.array(resol_list)[wsel]
            return s_list, z_list, logN_list, resol_list



        k_list = []
        id_list = []
        for i, s in enumerate(series_split):
            #print(s, 'systs_like')
            #series_o = list(set(series_split) - set([s]))
            trans = trans_parse(s)
            #z_int = np.arange(z_start, z_end, dz)
            if s in likes.keys():
                #print(likes[s])
                z_int = z_likes[s]
                #print(s)
                #plt.plot(z_int, likes[s])
                #w = np.where(likes[s]>sigma)
                #plt.plot(z_int[w], likes[s][w])
                #print(len(w[0]))
                """
                for s_o in series_o:
                    trans_o = trans_parse(s_o)
                    for t in trans_o:
                        print(t)
                        t_o = np.interp(x_w, spec._t['x'], spec._t[t])
                        print(likes_o)
                """
                #p0, _ = find_peaks(likes[s][w], distance=distance)
                #plt.scatter(z_int[w][p0], likes[s][w][p0])

                p0, _ = find_peaks(likes[s], distance=distance)
                pw = np.where(likes[s][p0]>sigma)
                p0 = p0[pw]
                #plt.scatter(z_int[p0], likes[s][p0])

                """
                # Check if likelihood peaks are higher than those of all other
                # transitions at those wavelengths
                x_w = np.array([to_x(z_int[w][p0], t) for t in trans])
                #for x_wi in x_w:
                #    plt.scatter(x_wi, likes[s][w][p0])
                t_all = np.array([])
                #for so in np.array(series_split):
                for so in likes.keys():
                    t_all = np.append(t_all, trans_parse(so))
                #t_all = np.ravel([trans_parse(so) for so in np.array(series_split)])
                liket = np.array([])
                for to in t_all:
                    liket = np.append(liket, np.interp(x_w, spec._t['x'], spec._t[to]))
                    #print(to, liket)
                liket = np.reshape(liket, (len(t_all), len(np.ravel(x_w))))
                liket[np.isnan(liket)] = -9999
                check = [t in trans for t in t_all[np.argmax(liket, axis=0)]]
                sel = np.prod(np.reshape(check, (len(trans), len(p0))), axis=0)
                #print(sel)
                #print(p0)
                #print(p0[np.where(sel)])
                """
                #print(z_int[w][p0])
                #x_w = to_x(z_int[w][p0], trans[0])
                x_w = to_x(z_int[p0], trans[0])
                #print(x_w)
                psel = []
                for p, x in zip(p0, x_w):
                    c = np.abs(x - spec._t['x']).argmin()
                    if spec._t['like_'+s][c]==spec._t['like'][c]:
                        psel.append(p)


                p = psel

                for ssub in s.split(':'):
                    if k_list == []:
                        s_list = [ssub]*len(p)

                        #z_list = z_int[w][p]
                        z_list = z_int[p]
                        logN_list = [logN]*len(p)
                        resol_list = [resol]*len(p)

                        # Check that components do not fall in masked regions
                        """
                        sel = []
                        for ssel,zsel in zip(s_list, z_list):
                            xint = []
                            for t in trans_parse(ssel):
                                xsel = to_x(zsel, t)
                                xint.append(np.interp(xsel, spec._t['x'],
                                                      spec._t['y']))
                            sel.append(not np.any(np.isnan(xint)))
                        wsel = np.where(sel)[0]
                        s_list = np.array(s_list)[wsel]
                        z_list = np.array(z_list)[wsel]
                        logN_list = np.array(logN_list)[wsel]
                        resol_list = np.array(resol_list)[wsel]
                        """
                        s_list, z_list, logN_list, resol_list = \
                            _nan_check(spec, s_list, z_list, logN_list, resol_list)

                        if len(s_list)>0:
                            self._systs_prepare(append)
                            id_list = self._systs_add(s_list, z_list, logN_list,
                                                      resol_list=resol_list)
                            self._spec_update()
                        else:
                            id_list = []
                    else:
                        s_list = [ssub]*len(z_list)
                        logN_list = [logN]*len(z_list)
                        resol_list = [resol]*len(z_list)

                        s_list, z_list, logN_list, resol_list = \
                            _nan_check(spec, s_list, z_list, logN_list, resol_list)

                        self._systs_prepare(append)
                        id_list = self._systs_add(s_list, z_list, logN_list,
                                                  resol_list=resol_list, k_list=k_list)
                        self._mods_recreate()
                        self._spec_update()
            else:
                id_list = []
            if i == 0:
                k_list = ['lines_voigt_%i_z' % id for id in id_list]
                #print(k_list)
        #plt.show()
        if compressed:
            systs._compress()
        #print('end')
        return 0


    def systs_merge(self, to_row=0, from_rows=[1]):
        """ @brief Merge a system into the current system
        @details Merged systems appear as a single entry in the compressed
        system table.
        @param num1 Row of the current system
        @param num2 Row of the system to be merged
        @return 0
        """

        try:
            to_row = int(to_row)-1
            from_rows = np.array(from_rows)-1
        except:
            logging.error(msg_param_fail)
            return 0

        t = self.sess.systs.t
        #print(t['z'][from_rows])
        z_app = np.append(t['z'][from_rows], t['z'][to_row])
        logN_app = np.append(t['logN'][from_rows], t['logN'][to_row])

        t['z'][to_row] = np.average(z_app, weights=10**logN_app)
        t['logN'][to_row] = np.log10(np.sum(10**logN_app))
        t.remove_rows(from_rows)
        t.sort(['z','id'])

        return 0

    def systs_new_from_lines(self, series='Ly-a', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=np.inf,
                             refit_n=0, chi2rav_thres=1e-2, max_nfev=max_nfev_def,
                             append=True):
        """ @brief New systems from line list
        @details Add and fit Voigt models to a line list, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param refit_n Number of refit cycles
        @param chi2rav_thres Average chi2r variation threshold between cycles
        @param max_nfev Maximum number of function evaluation
        @param append Append systems to existing system list
        @return 0
        """

        try:
            #series = series.replace(';',',')
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            resol = None if resol in [None, 'None'] else float(resol)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            self._refit_n = int(refit_n)
            self._chi2rav_thres = float(chi2rav_thres)
            self._max_nfev = float(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        check, resol = resol_check(self.sess.spec, resol)
        if not check: return 0

        systs = self.sess.systs

        recompress = False
        compressed = False
        if systs is not None and systs._compressed:
            recompress = True
            systs._compress()


        for s in series.split(';'):
            z_list, y_list = self._lines_cands_find(s, z_start, z_end, dz)
            z_list, logN_list, _ = self.sess.lines._cands_find2(s, z_start, z_end,
                                                            dz, logN=logN is None)

            if len(z_list) == 0:
                logging.warning("I've found no candidates!")
                return 0

            s_list = [s]*len(z_list)
            resol_list = [resol]*len(z_list)


            self._systs_prepare(append)
            #self._logN_guess(series, z_list[0], b, resol)
            #logN_list = self._systs_guess(series_list, z_list)
            self._systs_add(s_list, z_list, logN_list, resol_list=resol_list)
            #self._systs_fit()
            self._systs_cycle()
            self._spec_update()

        if compressed:
            systs._compress()

        return 0

    def syst_new_from_resids_new(self, series='Ly-a', z_start=0, z_end=6,
                             dz=1e-4, logN=logN_def, b=b_def, resol=resol_def,
                             chi2r_thres=np.inf, dlogN_thres=0.5,
                             max_nfev=max_nfev_def, append=True):
        """ @brief New system from residuals
        @details Add and fit a Voigt model from the strongest residual of
        previously fitted models in the neighborhood.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param logN Guess (logarithmic) column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param dlogN_thres Column density error threshold to accept the fitted model
        @param max_nfev Maximum number of function evaluation
        @param append Append systems to existing system list
        @return 0
        """

        try:
            series = series.replace(';',',')
            z_start = float(z_start)
            z_end = float(z_end)
            if series == 'unknown':
                z_start = 0
                z_end = np.inf
            dz = float(dz)
            if logN is not None:
                logN = float(logN)
            b = float(b)
            self._chi2r_thres = float(chi2r_thres)
            self._dlogN_thres = float(dlogN_thres)
            resol = float(resol)
            max_nfev = int(max_nfev)
            append = str(append) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        self.gauss_convolve(std=4, input_col='deabs', output_col='deabs_conv')
        sess = self.peaks_find(col='deabs_conv', kind='min', kappa=3.0, new_sess=True)


        #z_cand =

        return sess

    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def _systs_new_from_erf(self, series='Ly-a', col='y', z_start=0, z_end=6,
                            sigma=1, distance=3, append=True):
        """ @brief New systems from error function
        @details TBD
        @param series Series of transitions
        @param col Column to apply the likelihood
        @param z_start Start redshift
        @param z_end End redshift
        @param sigma Significance of absorbers (in units of the local error)
        @param distance Distance between systems in pixels
        @param append Append systems to existing system list
        @return 0
        """

        try:
            z_start = float(z_start)
            z_end = float(z_end)
            sigma = float(sigma)
            distance = None if distance in [None, 'None'] else float(distance)
        except:
            logging.error(msg_param_fail)
            return 0

        modul = 10
        distance = 3
        self.systs_new_from_like(series=series, col=col, z_start=z_start,
                                 z_end=z_end, modul=modul, sigma=sigma,
                                 distance=distance, append=append)

        return 0

    def _series_fit(self, series, zem, z_start=None, z_end=None, sigma=2, iter_n=3):

        def z_check(zem, z_start, z_end, s):
            if zem != None:
                if 'Ly_a' not in trans_parse(s):
                    z_start = (1+zem)*xem_d['Ly_a']/xem_d[trans_parse(s)[-1]]-1
                else:
                    z_start = (1+zem)*xem_d['Ly_b']/xem_d['Ly_a']-1
                z_end = zem
            return z_start, z_end

        def resid_peaks(std=10):
            spec = self.sess.spec
            #spec._gauss_convolve(std=std, input_col='deabs', output_col='deabs_mod')
            spec.t['deabs_mod'] = np.maximum(1, spec.t['deabs'])
            p, _ = find_peaks(spec.t['deabs_mod'], prominence=0.05)
            #print(np.array(spec._t['x'][p]))

        for s in series.split(';'):
            z_start, z_end = z_check(zem, z_start, z_end, s)
            self._systs_new_from_erf(series=s, z_start=z_start, z_end=z_end,
                                     sigma=sigma)
        self.systs_fit(refit_n=1)#, max_nfev=0)
        for i in range(iter_n):
            self.sess.systs._freeze_pars()
            for s in series.split(';'):
                z_start, z_end = z_check(zem, z_start, z_end, s)
                self._systs_new_from_erf(series=s, col='deabs',
                                         z_start=z_start, z_end=z_end, sigma=sigma)
            self.systs_fit(refit_n=0)#, max_nfev=0)
            self.sess.systs._unfreeze_pars()
            resid_peaks()
        if iter_n > 0:
            self.systs_fit(refit_n=1)
        #plt.show()


    def lya_fit(self, zem=None, z_start=None, z_end=None, sigma=1, iter_n=3):
        """ @brief Fit the Lyman-alpha forest
        @details The recipe identifies Lyman-alpha absorbers using the
        likelihood method and fits them. The procedure is iterated on residuals
        of the fit to improve it.
        @param zem Emission redshift
        @param z_start Start redshift (ignored if zem is specified)
        @param z_end End redshift (ignored if zem is specified)
        @param sigma Significance of absorbers (in units of the local error)
        @param iter_n Number of iterations on residuals
        @return 0
        """

        try:
            zem = None if zem in [None, 'None'] else float(zem)
            z_start = None if z_start in [None, 'None'] else float(z_start)
            z_end = None if z_end in [None, 'None'] else float(z_end)
            sigma = float(sigma)
            iter_n = int(iter_n)
        except:
            logging.error(msg_param_fail)
            return 0

        if zem != None:
            z_start = (1+zem)*xem_d['Ly_b']/xem_d['Ly_a']-1
            z_end = zem

        self._series_fit('Ly_a', zem, z_start, z_end, sigma, iter_n)
        #plt.show()
        return 0


    #@arg_fix(arg_mapping={'thres': 'sigma'})
    def red_fit(self,
                series='CIV;SiIV:CIV;SiII_1526:CIV;AlIII;MgII_2796,MgII_2803;'\
                       +'FeII_2382,FeII_2600:MgII_2796,MgII_2803',
                zem=None, z_start=None, z_end=None, sigma=1, iter_n=3):
        """ @brief Fit the red part of the spectrum forest
        @details The recipe identifies Lyman-alpha absorbers using the
        likelihood method and fits them. The procedure is iterated on residuals
        of the fit to improve it.
        @param zem Emission redshift
        @param z_start Start redshift (ignored if zem is specified)
        @param z_end End redshift (ignored if zem is specified)
        @param sigma Significance of absorbers (in units of the local error)
        @param iter_n Number of iterations on residuals
        @return 0
        """

        try:
            zem = None if zem in [None, 'None'] else float(zem)
            z_start = None if z_start in [None, 'None'] else float(z_start)
            z_end = None if z_end in [None, 'None'] else float(z_end)
            sigma = float(sigma)
            iter_n = int(iter_n)
        except:
            logging.error(msg_param_fail)
            return 0

        """
        for s in ['CIV','MgII_2796,MgII_2803','SiIV']:
            t = trans_parse(s)
            if zem != None:
                z_start = (1+zem)*xem_d['Ly_a']/xem_d[t[-1]]-1
                z_end = zem
        """
        self._series_fit(series, zem, z_start, z_end, sigma, iter_n)
            #self._systs_new_from_erf(series=s, z_start=z_start, z_end=z_end,
            #                         sigma=sigma)

        return 0
