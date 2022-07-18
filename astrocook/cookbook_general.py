from .functions import _gauss, expr_eval, running_mean, running_rms
from .message import *
from .vars import *
import ast
from astropy import table as at
from copy import deepcopy as dc
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys
from tqdm import tqdm

class CookbookGeneral(object):
    """ Cookbook of general utilities
    @details This cookbook contains utilities to manipulate sessions, mask the
    spectrum, estimate spectral quality parameters, and perform basic operations
    like rebinning and convolution.
    """

    def __init__(self):
        super(CookbookGeneral, self).__init__()

    def bin_zap(self, x):
        self.sess.spec._zap(xmin=x, xmax=None)




    def combine(self, name='*_combined', _sel=''):
        """ @brief Combine two or more sessions
        @details Create a new session combining the spectra from two or more
        other sessions.

        The recipe collects all the bins from the original spectra and puts them
        all together in the new spectrum. The bins retain their original size
        (defined by `xmin` and `xmax`), so they may overlap in the final
        spectrum. By default, they are ordered by ascending `x`.

        All other structures from the original sessions (line lists, etc.) are
        not propagated to the new one.

        N.B. To select sessions, either click on the session window or provide
        a list through the hidden parameter `_sel`.
        @param name Name of the output session
        @return Combined session
        """
        name_in = name
        #sel = self._tab._get_selected_items()
        sel = self.sess._gui._sess_item_sel
        sess_list = self.sess._gui._sess_list

        """
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            try:
                sel = [int(s) \
                       for s in _sel.replace('[','').replace(']','').split(',')]
            except:
                pass
        if sel == []:
            sel = range(len(sess_list))
        self._gui._sess_item_sel = sel
        """
        sel = _sel

        struct_out = {}
        for struct in sess_list[sel[0]].seq:
            struct_out[struct] = dc(getattr(sess_list[sel[0]], struct))


        if name_in[0] == '*':
            name = sess_list[sel[0]].name

        logging.info("Combining sessions %s..." % ', '.join(str(s) for s in sel))
        for s in sel[1:]:
            #spec._t = at.vstack([spec._t, self._gui._sess_list[s].spec._t])

            for struct in sess_list[s].seq:
                if getattr(sess_list[s], struct) != None:
                    if struct_out[struct] != None:
                        struct_out[struct]._append(
                            getattr(sess_list[s], struct))
                    else:
                        struct_out[struct] = dc(getattr(sess_list[s], struct))

            if name_in[0] == '*':
                name += '_' + sess_list[s].name

        struct_out['spec']._t.sort('x')
        if name_in[0] == '*':
            name += name_in[1:]
        from .session import Session
        sess = Session(gui=self.sess._gui, name=name, spec=struct_out['spec'],
                       nodes=struct_out['nodes'], lines=struct_out['lines'],
                       systs=struct_out['systs'])
        return sess


    def equalize(self, xmin, xmax, _sel='', cont=True):
        """ @brief Equalize two sessions
        @details Equalize the spectrum of two sessions, based on their flux
        ratio within a wavelength window.

        By default, the last-selected spectrum is equalized to the
        first-selected one (which is left unchanged). Equalization is done in
        place, without creating a new session.

        To compute the rescaling factor, the recipe takes the medians of the `y`
        columns of the two spectra between `xmin` and `xmax`. The `y` and `dy`
        columns of the second spectrum are then multiplied by $$
        \\textrm{med}($$`y`$$_1)/\\textrm{med}($$`y`$$_2)$$.

        N.B. To select sessions, either click on the session window or provide
        a list through the hidden parameter `_sel`.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return 0
        """

        try:
            xmin = float(xmin) * au.nm
            xmax = float(xmax) * au.nm
        except ValueError:
            logging.error(msg_param_fail)
            return None

        """
        sel = self._gui._sess_item_sel
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            sel = [int(s) \
                for s in _sel.replace('[','').replace(']','').split(',')]
        self._gui._sess_item_sel = sel
        """
        sel = _sel
        logging.info("Equalizing session %i to session %i... "
                     % (sel[1], sel[0]))

        for i,s in enumerate(sel):
            sess = self.sess._gui._sess_list[s]
            w = np.where(np.logical_and(sess.spec.x>xmin, sess.spec.x<xmax))[0]
            if len(w)==0:
                logging.error("I can't use this wavelength range for "
                              "equalization. Please choose a range covered by "
                              "both sessions.")
                return(0)
            if i == 0:
                f = np.nanmedian(sess.spec.y[w]).value
                #print(np.median(sess.spec.y[w]))
            else:
                f = f/np.nanmedian(sess.spec.y[w]).value
                #print(np.median(sess.spec.y[w]), f)
                logging.info("Equalization factor: %3.4f." % f)
                sess.spec.y = f*sess.spec.y
                sess.spec.dy = f*sess.spec.dy
                if cont and 'cont' in sess.spec._t.colnames:
                    sess.spec._t['cont'] = f*sess.spec._t['cont']

        return 0


    def dx_est(self):
        """ @brief Estimate bin size in x
        @details Compute statistics on xmax-xmin, to determine the typical
        binsize in wavelength and velocity space.
        @return 0
        """

        spec = self.sess.spec
        dx = spec.t['xmax'][1:-1]-spec.t['xmin'][1:-1]
        dx_mean = np.mean(dx)
        dx_std = np.std(dx)
        logging.info('Distribution of dx: %.6f±%.6f' % (dx_mean, dx_std))

        return 0


    def feature_zap(self, xmin, xmax):
        self.sess.spec._zap(xmin, xmax)


    def flux_ccf(self, col1='y', col2='y', dcol1='dy', dcol2='dy', vstart=-20,
                 vend=20, dv=1e-1):
        """@brief Compute the CCF
        @details Convolve the cross-correlation function (CCF) between two
        spectrum columns.

        The recipe is designed to work on flux densities. Typically, the first
        column is `y` and the second column contains the flux density from a
        different spectrum with the same wavelength binning. The second columns
        can also be `y`: in this case, the recipe computes the auto-correlation
        instead of the cross-correlation.

        The CCF is computed in velocity space, shifting `col2` with respect to
        `col1` within the range `vstart`-`vend` and with step @dv. The columns
        are resampled while shifting, to accomodate for values of @dv much
        smaller than the spectrum bin size.

        The CCF is saved in a NumPy binary file `SESS_ccf.npy`, with `SESS` the
        name of the session.
        @param col1 First column
        @param col2 Second column
        @param dcol1 Error for first column
        @param dcol2 Error for second column
        @param vstart Start velocity
        @param vend End velocity
        @param dv Velocity step
        @return 0
        """

        try:
            vstart = float(vstart) * au.km/au.s
            vend = float(vend) * au.km/au.s
            dv = float(dv) * au.km/au.s
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        v_shift, ccf = self.sess.spec._flux_ccf(col1, col2, dcol1, dcol2, vstart,
                                                vend, dv)
        logging.info("CCF statistics: minimum %3.4f, maximum %3.4f, mean %3.4f." \
                     % (np.min(ccf), np.max(ccf), np.mean(ccf)))
        with open(self.sess.name+'_ccf.npy', 'wb') as f:
            np.save(f, v_shift)
            np.save(f, ccf)

        return 0


    def flux_ccf_stats(self, n=1e1, col1='y', col2='y', dcol1='dy', dcol2='dy',
                       vstart=-20, vend=20, dv=1e-1, fit_hw=1.):
        """@brief Compute statistics of the CCF
        @details Compute statistics for the peak of the cross-correlation
        function (CCF) by bootstrapping a number of realizations for the
        spectrum.

        Realizations are created by selecting entries at random, preserving
        wavelength order and rejecting duplicates (compare with Peterson et al.
        1998).

        The recipe computes the CCF between the original flux and the flux of
        each realization. A gaussian is fit to the CCF within a window around
        0 (in velocity space) to determine the position of the peak. The
        distribution of peak positions is saved in a NumPy binary file
        `SESS_ccf_stats.npy`, with `SESS` the name of the session.
        @param n Number of realizations
        @param col1 First column
        @param col2 Second column
        @param dcol1 Error for first column
        @param dcol2 Error for second column
        @param vstart Start velocity (km/s)
        @param vend End velocity (km/s)
        @param dv Velocity step (km/s)
        @param fit_hw Half-window used for fitting the CCF (km/s)
        @return 0
        """

        try:
            n = int(n)
            vstart = float(vstart) * au.km/au.s
            vend = float(vend) * au.km/au.s
            dv = float(dv) * au.km/au.s
            fit_hw = float(fit_hw) * au.km/au.s
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        peaks = np.array([])
        shifts = np.array([])
        #for i in range(n):
        for i in enum_tqdm(range(n), len(range(n)), "cookbook_general: Bootstrapping"):
            spec = dc(self.sess.spec)
            rng = np.random.default_rng()
            rint = rng.integers(0, len(spec.t), size=len(spec.t))
            sel = np.sort(np.unique(rint))
            #sel = np.unique(rint)
            spec._t = spec._t[sel]
            #plt.plot(spec._t['x'], spec._t['y'])
            v_shift, ccf = spec._flux_ccf(col1, col2, dcol1, dcol2, vstart,
                                          vend, dv)


            v_shiftmax = v_shift[np.argmax(ccf)]
            try:
            #    ciao
            #except:
                p0 = [1., v_shiftmax, 1.]
                fit_sel = np.logical_and(v_shift>v_shiftmax-fit_hw.value,
                                         v_shift<v_shiftmax+fit_hw.value)
                #plt.plot(v_shift[fit_sel], ccf[fit_sel], linestyle=':')
                coeff, var_matrix = curve_fit(_gauss, v_shift[fit_sel], ccf[fit_sel], p0=p0)
                fit = _gauss(v_shift[fit_sel], *coeff)
                #perr = np.sqrt(np.diag(var_matrix))
                peak, shift = coeff[:2]
                #print(peak, shift)
            except:
                peak, shift = np.nan, np.nan
            peaks = np.append(peaks, peak)
            shifts = np.append(shifts, shift)
        #print(peaks, shifts)
            #logging.info("CCF statistics: minimum %3.4f, maximum %3.4f, "
            #             "mean %3.4f, shift %3.4f." \
            #             % (np.min(ccf), np.max(ccf), np.mean(ccf), shift))

        logging.info("Peak: %3.4f±%3.4f; shift: %3.4f±%3.4f" \
            % (np.nanmean(peaks), np.nanstd(peaks), np.nanmean(shifts), np.nanstd(shifts)))
        #plt.show()
        with open(self.sess.name+'_ccf_stats.npy', 'wb') as f:
            np.save(f, peaks)
            np.save(f, shifts)
        #plt.show()
        return 0


    def gauss_convolve(self, std=20.0, input_col='y', output_col='conv'):
        """@brief Convolve with gaussian
        @details Convolve a spectrum column with a gaussian profile.

        The convolution is computed in velocity space, using the Fast Fourier
        Transform.
        @param std Standard deviation of the gaussian (km/s)
        @param input_col Input column
        @param output_col Output column
        @return 0
        """

        try:
            std = float(std) * au.km/au.s
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._gauss_convolve(std, input_col, output_col)
        return 0

    def mask(self, col='mask', cond='', new_sess=True, masked_col='x'):
        """ @brief Mask the spectrum
        @details Create a mask applying a specified condition to the spectrum
        bins.

        The expression in `cond` is translated into a boolean condition by the
        [`ast`](https://docs.python.org/3/library/ast.html) module. Expressions
        like `c>10` or `1500<c<2000` are supported, where `c` is a column of the
        spectrum.

        The condition is checked on all spectrum bins and a new column `col` is
        populated with the results. No information is deleted from the input
        spectrum. If `new_sess` is `True`, a new session is created, containing
        a masked version of the input spectrum. In this masked spectrum, the
        column `y`, `dy`, and optionally `cont` are set to `numpy.nan` in all
        bins where the condition is false. All other structures from the
        original session (line lists, etc.) are not propagated to the new one.
        @param col Column with the mask
        @param cond Condition
        @param new_sess Create a new session from masked spectrum
        @return Session with masked spectrum
        """

        spec = self.sess.spec

        for c in sorted(spec._t.colnames, key=len, reverse=True):
            cond = cond.replace(c, str(list(np.array(spec._t[c]))))
            #print(c, cond)
        mask = expr_eval(ast.parse(cond, mode='eval').body)

        if col not in spec._t.colnames:
            logging.info("I'm adding column %s." % col)
        else:
            logging.info("I'm updating column %s." % col)
        spec._t[col] = mask

        if new_sess:
            spec_out = dc(spec)
            #spec_out._t['x'][~mask] = np.nan
            for c in ['y', 'dy', 'cont']:
                if c in spec_out._t.colnames:
                    spec_out._t[c] = spec_out._t[c].astype(float)
                    spec_out._t[c][~mask] = np.nan #
            #spec_out._t = .spec._t[mask]
            from .session import Session
            new = Session(gui=self.sess._gui, name=self.sess.name+'_'+col,
                          spec=spec_out)
            return new

        return 0

    def part_extract(self, zem, part='blue'):
        """ @brief Extract blue or red part
        @details Extract blue or red part, based on the emission redshift.

        The recipe computes the observed wavelength of the Lyman alpha emission
        line as $$(1+$$`zwm`$$)\times 121.567\textrm{ nm}$$ and extracts the
        region bluewards or redwards of this wavelength into a new session,
        based on the value of `part`.

        @param zem Emission redshift
        @param part Either `blue` or `red`
        """

        try:
            zem = float(zem)
        except ValueError:
            logging.error(msg_param_fail)
            return None

        if part not in ['blue', 'red']:
            logging.error(msg_param_fail)
            return None

        if part == 'blue':
            return self.region_extract(0, (1+zem)*121.567)
        else:
            return self.region_extract((1+zem)*121.567, np.infty)


    def rebin(self, xstart=None, xend=None, dx=10.0, xunit=au.km/au.s,
              norm=False, filling=np.nan):
        """ @brief Re-bin spectrum
        @details Apply a new binning to a spectrum, with a constant bin size.

        The algorithm for re-binning is described in
        [Cupani et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016SPIE.9913E..1TC/abstract).
        It properly weights the flux contributions from the old bins to the
        new ones, also when the former overlap with each other (as it happens
        when several exposures of the same object are combined into a single
        spectrum).

        The new grid is designed to fully cover the original range of
        the spectrum (when `xstart` and `xend` are `None`) or a specified range
        (useful when different spectra must be re-binned to the same grid). It
        is defined in either wavelength or velocity space, as specified by the
        chosen `xunit`. Any gap in the original binning are filled with a
        specified `filling` value, to ensure that the final grid is equally
        spaced.

        Columns `y` and `dy` of the input spectrum are both re-binned to the new
        grid. If a column `cont` is present and `norm` is `True`, `y` and `dy`
        are normalized to `cont` in the re-binned spectrum.

        A new session is created with the re-binned spectrum. All other
        structures from the original session (line lists, etc.) are not
        propagated to the new one.
        @param xstart Start wavelength (nm)
        @param xend End wavelength (nm)
        @param dx Step in x
        @param xunit Unit of wavelength or velocity
        @param norm Return normalized spectrum, if continuum exists
        @param filling Value to fill region without data
        @return Session with rebinned spectrum
        """

        try:
            xstart = None if xstart in [None, 'None'] else float(xstart)
            xend = None if xend in [None, 'None'] else float(xend)
            dx = float(dx)
            xunit = au.Unit(xunit)
            norm = str(norm) == 'True'
            filling = float(filling)
        except ValueError:
            logging.error(msg_param_fail)
            return None

        """
        sel = self.sess._gui._sess_item_sel
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            sel = [int(s) \
                for s in _sel.replace('[','').replace(']','').split(',')]
        if sel == [] or len(sel)>1:
            sel = [self.sess._gui._panel_sess._tab.GetItemCount()-1]
        self.sess._gui._sess_item_sel = sel
        """

        #print(self.sess)
        # A deep copy is created, so the original spectrum is preserved
        spec_in = dc(self.sess.spec)

        cont = 'cont' in spec_in.t.colnames
        if not norm or (norm and not cont):
            if not cont:
                logging.warning("I can't find continuum to normalize the "
                                "spectrum. Using non-normalized y column "
                                "instead.")
            y = spec_in.y
            dy = spec_in.dy
        else:
            y = spec_in.y/spec_in.t['cont']
            dy = spec_in.dy/spec_in.t['cont']

        spec_out = spec_in._rebin(xstart, xend, dx, xunit, y, dy, filling)
        if cont:
            if not norm:
                spec_out.t['cont'] = np.interp(
                    spec_out.x.to(au.nm).value,spec_in.x.to(au.nm).value,
                    spec_in.t['cont']) * spec_out.y.unit
            else:
                spec_out.t['cont'] = np.ones(len(spec_out.t)) * spec_out.y.unit

        # Create a new session
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_rebinned',
                      spec=spec_out)
        return new


    def region_extract(self, xmin, xmax, verbose=True):
        """ @brief Extract region
        @details Extract a spectral region and create a new session from it.

        The region includes not only the chunk of spectrum between `xmin` and
        `xmax`, but also all the lines and the absorption systems in the same
        range, which are propagated to the new session.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return Spectral region
        """
        try:
            xmin = float(xmin) * au.nm
            xmax = float(xmax) * au.nm
        except ValueError:
            logging.error(msg_param_fail)
            return None

        if xmin > xmax:
            temp = xmin
            xmin = xmax
            xmax = temp
            logging.warning(msg_param_swap)

        kwargs = {'path': self.sess.path, 'name': self.sess.name}
        for s in self.sess.seq:
            try:
                struct = getattr(self.sess, s)._region_extract(xmin, xmax)
                if struct is None:
                    logging.warning(msg_empty(s))
                else:
                    kwargs[s] = struct
            except:
                logging.debug("Attribute %s does not support region "
                              "extraction." % s)
                try:
                    kwargs[s] = getattr(self.sess, s)
                except:
                    logging.debug(msg_attr_miss(s))
                    kwargs[s] = None
        if kwargs['spec'] != None:
            from .session import Session
            new = Session(**kwargs)
            new._gui = self.sess._gui
        else:
            new = None
        if 'systs' in self.sess.seq and self.sess.systs != None \
            and new.systs != None:

            # This is needed instead of a simple deepcopy because
            # GUIPanelSession does not support pickling
            #old = dc(self.sess)
            old = Session(self.sess._gui)
            for d in self.sess.__dict__:
                if d != '_gui' and d != 'cb' and d != 'log' and d != 'defs':
                    old.__dict__[d] = dc(self.sess.__dict__[d])
                if d == 'defs':
                    setattr(old, d, getattr(self.sess, d))
                    for dd in getattr(self.sess, d).__dict__:
                        if dd == '_gui':
                            getattr(old, d).__dict__[dd] = self.sess._gui
                        else:
                            getattr(old, d).__dict__[dd] \
                                = dc(getattr(self.sess, d).__dict__[dd])
                    #print(getattr(self.sess, d).__dict__)
            old.__dict__['cb'] = self.sess.__dict__['cb']

            self.sess = new
            self._mods_recreate()
            self.sess = old

        return new


    def resol_est(self, px=3, update=True):
        """ @brief Estimate resolution
        @details Assign a resolution to spectral bins, assuming that the
        spectrum is designed to have a fixed number of bins per resolution
        element.

        This recipe is useful to populate the `resol` column in a spectrum
        (needed to fit the absorption systems) when it is empty, and information
        about the original sampling of the data is available. It does *not* try
        to infer the resolution from, e.g., the width of unresolved spectral
        feature.
        @param px Number of bins per resolution element
        @param update Update column 'resol'
        @return 0
        """

        try:
            px = int(px)
            update = str(update) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._resol_est(px, update)

        return 0

    def _rm_est(self, hwindow=100):
        """ @brief Compute running mean
        @details Compute running mean of the flux.
        @param hwindow Half-window size in pixels for running mean
        @return 0
        """

        try:
            hwindow = int(hwindow)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec

        y_rm = running_mean(spec._t['y'], h=hwindow)
        if 'y_rm' not in spec._t.colnames:
            logging.info("I'm adding column 'y_rm'.")
        else:
            logging.warning("I'm updating column 'y_rm'.")
        spec._t['y_rm'] = at.Column(y_rm, dtype=float)

        return 0

    def _rfz_set(self, z=0.0):
        """ @brief Set redshift for a rest-frame spectrum
        @details Set a redshift value for a spectrum shifted to rest frame.
        @param z redhisft
        @return 0
        """

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._rfz_man = z
            except:
                logging.debug(msg_attr_miss(s))

        return 0


    def rms_est(self, hwindow=100):
        """ @brief Estimate error from RMS
        @details Estimate flux error by computing the root-mean-square (RMS) of
        the flux within a running window.

        The RMS is computed over `y` values and is saved in `y_rms`. It may be
        useful to compare the latter with `dy` to check that the formal error is
        consistent with the actual dispersion of `y` values.
        @param hwindow Half-size in pixels of the running window
        @return 0
        """

        try:
            hwindow = int(hwindow)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec

        y_rm = running_mean(spec._t['y'], h=5)
        y_rms = running_rms(spec._t['y'], y_rm, h=hwindow)
        if 'y_rms' not in spec._t.colnames:
            logging.info("I'm adding column 'y_rms'.")
        else:
            logging.warning("I'm updating column 'y_rms'.")
        spec._t['y_rms'] = at.Column(y_rms, dtype=float)

        return 0


    def snr_est(self):
        """ @brief Estimate the SNR
        @details Estimate the signal-to-noise ratio per pixel.

        A `snr` column is populated with `y`/`dy` ratios computed for all
        spectrum bins.
        @return 0
        """

        spec = self.sess.spec
        if 'snr' not in spec._t.colnames:
            logging.info("I'm adding column 'snr'.")
        else:
            logging.info("I'm updating column 'snr'.")

        spec._t['snr'] = spec.y/spec.dy

        return 0


    def telluric_mask(self, shift=0, thres=0.99, apply=True):
        """ @brief Mask telluric absorption
        @details Mask spectral regions affected by telluric absorptions.

        The regions were determined by Tobias M. Schmidt from ESPRESSO data and
        are saved in `telluric.dat`. They are resampled into the current `x`
        grid and used to populate a `telluric` column, which is set to `1`
        inside the regions and to `0` elsewhere.

        If `apply` is `True`, `y` is set to `numpy.nan` in all bins where
        `telluric` is 1.
        @param shift Shift to the heliocentric frame (km/s)
        @param thres Threshold to cut telluric lines in the model (normalized to
        continuum)
        @param apply Apply mask to flux
        @return 0
        """

        try:
            shift = float(shift)
            thres = float(thres)
            apply = str(apply) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        if 'telluric' not in spec._t.colnames:
            logging.info("I'm adding column 'telluric'.")
        else:
            logging.info("I'm updating column 'telluric'.")

        #p = '/'.join(pathlib.PurePath(os.path.realpath(__file__)).parts[0:-1]) + '/../'
        #telluric = ascii.read(pathlib.Path(p+'/telluric.dat'))
        #telluric = fits.open(pathlib.Path(p+'/telluric.fits'))[1].data
        x = np.array(telluric['lam'], dtype=float) * (1+shift/aconst.c.to(au.km/au.s).value)
        mask = np.array([t<thres for t in telluric['trans_ma']], dtype=float)
        tell = np.interp(spec._t['x'].to(au.nm).value, x, mask)
        spec._t['telluric'] = np.array(tell!=1, dtype=bool)

        if apply:
            spec._t['y'][np.where(np.array(tell==1, dtype=bool))] = np.nan

        return 0
