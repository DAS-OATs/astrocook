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
        @details Combine the spectra from two or more sessions. The sessions are
        either selected on the Sessions window or provided as a list through the
        @_sel parameter. A new session is created, with a new spectrum
        containing all entries from the spectra of the combined sessions. Other
        objects from the sessions (line lists, etc.) are discarded.
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


    def deredden(self, ebv=0.03, rv=3.1):
        """@brief Deredden spectrum
        @details Deredden the spectrum using the parametrization by Cardelli,
        Clayton, and Mathis (1989) and O'Donnell (1994).
        @param ebv Color excess E(B-V)
        @param rv Ratio of total selective extinction R(V)=A(V)/E(B-V)
        @return 0
        """

        try:
            ebv = float(ebv)
            rv = float(rv)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._deredden(ebv, rv)
        return 0


    def equalize(self, xmin, xmax, _sel=''):
        """ @brief Equalize two sessions
        @details Equalize the flux level of one session to another one. The
        last-selected session is equalized to the first-selected one. The
        equalization factor is the ratio of the median flux within the
        specified wavelength interval.
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

        return 0


    def feature_zap(self, xmin, xmax):
        self.sess.spec._zap(xmin, xmax)


    def flux_ccf(self, col1='y', col2='y', dcol1='dy', dcol2='dy', vstart=-20,
                 vend=20, dv=1e-1):
        """@brief Compute the flux CCF
        @details Convolve the cross-correlation function between two columns
        with flux information. Typically, the second column contain the flux
        density from a different spectrum, but it can be equal to the first
        column (in which case auto-correlation is computed). Cross-correlation
        is computed in velocity space, within a given range and with a given
        step.
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
        #plt.show()

        return 0


    def flux_ccf_stats(self, n=1e1, col1='y', col2='y', dcol1='dy', dcol2='dy',
                       vstart=-20, vend=20, dv=1e-1, fit_hw=1.):
        """@brief Compute statistics of the flux CCF
        @details Compute statistics of the flux CCF by bootstrapping a number
        of realizations for the spectrum. Realizations are created by selecting
        entries at random, preserving wavelength order and rejecting duplicates.
        Compare with Peterson et al. 1998.
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
            v_shift, ccf = spec._flux_ccf(col1, col2, dcol1, dcol2, vstart,
                                          vend, dv)

            try:
            #    ciao
            #except:
                p0 = [1., 0., 1.]
                fit_sel = np.logical_and(v_shift>-fit_hw.value, v_shift<fit_hw.value)
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

            #logging.info("CCF statistics: minimum %3.4f, maximum %3.4f, "
            #             "mean %3.4f, shift %3.4f." \
            #             % (np.min(ccf), np.max(ccf), np.mean(ccf), shift))

        logging.info("Peak: %3.4f±%3.4f; shift: %3.4f±%3.4f" \
            % (np.nanmean(peaks), np.nanstd(peaks), np.nanmean(shifts), np.nanstd(shifts)))
        #plt.show()
        with open(self.sess.name+'_ccf_stats.npy', 'wb') as f:
            np.save(f, peaks)
            np.save(f, shifts)

        return 0


    def gauss_convolve(self, std=20.0, input_col='y', output_col='conv'):
        """@brief Convolve with gaussian
        @details Convolve a spectrum column with a gaussian profile using FFT
        transform.
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
        """ @brief Create a spectral mask
        @details Create a spectral mask by applying a given condition. The
        condition must be parsable by AST, with spectrum columns denoted by
        their names (e.g. 'x>400'). Optionally, a new session is created with
        the masked spectrum. Other objects from the old session (line lists,
        etc.) are discarded.
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
            #spec_out._t['x'][~mask] = np.nan #
            for c in ['y', 'dy', 'cont']:
                if c in spec_out._t.colnames: spec_out._t[c][~mask] = np.nan #
            #spec_out._t = .spec._t[mask]
            from .session import Session
            new = Session(gui=self.sess._gui, name=self.sess.name+'_'+col,
                          spec=spec_out)
            return new

        return 0


    def rebin(self, xstart=None, xend=None, dx=10.0, xunit=au.km/au.s,
              norm=True, filling=np.nan):
        """ @brief Rebin spectrum
        @details Rebin a spectrum with a given step. The step can be expressed
        in any unit of wavelength or velocity. Start and end wavelength may be
        specified, e.g. to align the rebinned spectrum to other spectra. If
        start or end wavelength are None, rebinning is performed from the first
        to the last wavelength of the input spectrum. A new session is created
        with the rebinned spectrum. Other objects from the old session (line
        lists, etc.) are discarded.
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
            spec_out.t['cont'] = 1

        # Create a new session
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_rebinned',
                      spec=spec_out)
        return new


    def region_extract(self, xmin, xmax):
        """ @brief Extract region
        @details The region between a minimum and a maximum wavelength is
        extracted from the data structures in the current session (these include
        the selected spectral range with all the lines and the absorption
        systems that fall within). A new session with the extracted data
        structures is created.
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
                kwargs[s] = getattr(self.sess, s)._region_extract(xmin, xmax)
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
        if 'systs' in self.sess.seq and self.sess.systs != None:

            # This is needed instead of a simple deepcopy because
            # GUIPanelSession does not support pickling
            #old = dc(self.sess)
            old = Session(self.sess._gui)
            for d in self.sess.__dict__:
                if d != '_gui' and d != 'cb' and d != 'log':
                    old.__dict__[d] = dc(self.sess.__dict__[d])
            old.__dict__['cb'] = self.sess.__dict__['cb']

            self.sess = new
            self._mods_recreate()
            self.sess = old

        return new


    def resol_est(self, px=3, update=True):
        """ @brief Estimate resolution
        @details Estimate spectral resolution assuming the spectrum has
        a fixed number of pixels per resolution element.
        @param px Number of pixels
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
        @details Estimate flux error by computing the running RMS of the flux.
        @param hwindow Half-window size in pixels for running mean
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
        @return 0
        """

        spec = self.sess.spec
        if 'snr' not in spec._t.colnames:
            logging.info("I'm adding column 'snr'.")
        else:
            logging.info("I'm updating column 'snr'.")

        spec._t['snr'] = spec.y/spec.dy

        return 0


    def telluric_mask(self, shift=0, apply=True):
        """ @brief Mask telluric absorption
        @details Mask telluric absorption
        @param shift Shift to the heliocentric frame (km/s)
        @param apply Apply mask to flux
        @return 0
        """

        try:
            shift = float(shift)
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
        x = np.array(telluric['WAVEL'], dtype=float) * (1+shift/aconst.c.to(au.km/au.s).value)
        mask = np.array([t!=0 for t in telluric['MASK']], dtype=float)
        tell = np.interp(spec._t['x'].to(au.nm).value, x, mask)
        spec._t['telluric'] = np.array(tell!=1, dtype=bool)

        if apply:
            spec._t['y'][np.where(np.array(tell==1, dtype=bool))] = np.nan

        return 0
