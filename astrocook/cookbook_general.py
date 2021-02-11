from .functions import expr_eval
from .message import *
from .vars import *
import ast
from copy import deepcopy as dc
import numpy as np
import sys
from tqdm import tqdm

class CookbookGeneral(object):
    """ Cookbook of general utilities
    """

    def __init__(self):
        super(CookbookGeneral, self).__init__()

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


    def feature_zap(self, xmin, xmax):
        self.sess.spec._zap(xmin, xmax)


    def flux_ccf(self, col1='y', col2='y', vstart=-20, vend=20, dv=1e-1):
        """@brief Compute the flux CCF
        @details Convolve the cross-correlation function between two columns
        with flux information. Typically, the second column contain the flux
        density from a different spectrum, but it can be equal to the first
        column (in which case auto-correlation is computed). Cross-correlation
        is computed in velocity space, within a given range and with a given
        step.
        @param y First column
        @param y Second column
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

        self.sess.spec._flux_ccf(col1, col2, vstart, vend, dv)
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

    def mask(self, col='mask', cond='', new_sess=True):
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
            spec_out._t = spec._t[mask]
            from .session import Session
            new = Session(gui=self.sess._gui, name=self.sess.name+'_'+col,
                          spec=spec_out)
            return new

        return 0


    def rebin(self, xstart=None, xend=None, dx=10.0, xunit=au.km/au.s,
              norm=True):
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
        @return Session with rebinned spectrum
        """

        try:
            xstart = None if xstart in [None, 'None'] else float(xstart)
            xend = None if xend in [None, 'None'] else float(xend)
            dx = float(dx)
            xunit = au.Unit(xunit)
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

        spec_out = spec_in._rebin(xstart, xend, dx, xunit, y, dy)
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
                if d != '_gui' and d != 'cb':
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

    def shift_from_rf(self, z=0):
        """ @brief Shift from rest frame
        @details Shift x axis from rest frame to the original frame.
        @param z Redshift to use for shifting
        @return 0
        """

        try:
            z = float(z)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        for s in self.sess.seq:
            try:
                z_to = z-getattr(self.sess, s)._rfz
                getattr(self.sess, s)._shift_rf(z_to)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def shift_to_rf(self, z=0):
        """ @brief Shift to rest frame
        @details Shift x axis to the rest frame.
        @param z Redshift to use for shifting
        @return 0
        """

        try:
            z = float(z)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._shift_rf(z)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def shift_bary(self, v=None):
        """ @brief Shift to barycentric frame
        @details Shift x axis to the barycentric frame of the solar system.
        @param v Velocity in the barycentric frame (km/s)
        @return 0
        """

        try:
            v = float(v)
        except:
            try:
                v = self.sess.spec.meta['v_bary']
            except ValueError:
                logging.error(msg_param_fail)
                return 0

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._shift_bary(v)
            except:
                logging.debug(msg_attr_miss(s))
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



    def x_convert(self, zem=0, xunit=au.km/au.s):
        """ @brief Convert x axis
        @details Convert the x axis to wavelength or velocity units. The x axis
        can be converted to any unit of wavelength or velocity (default: nm and
        km/s). The conversion applies to both the spectrum and the line list.
        When converting to and from velocity units, the zero point is set at
        (1+zem)λ_Lya (where λ_Lya = 121.567 nm is the rest-frame wavelength
        of the Lyman-alpha transition).
        @param zem Emission redshift
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        try:
            zem = float(zem)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        xunit = au.Unit(xunit)
        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._x_convert(zem, xunit)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def y_convert(self, yunit=au.electron/au.nm):
        """ @brief Convert y axis
        @details Convert the y axis to different units. The y axis can be
        expressed in different units depending on how it was calibrated
        (default: erg/(cm^2 s nm)). It can be converted to any unit of the same
        physical quantity. The conversion applies to both the spectrum and the
        line list.
        @param yunit Unit of flux density
        @return 0
        """

        yunit = au.Unit(yunit)

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._y_convert(yunit=yunit)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def y_scale(self, fact=1.0):
        """ @brief Scale y axis
        @details Scale the y axis by a constant factor. The spectrum and the
        line list are rescaled in place, without starting a new session.
        @param fact Multiplicative factor
        @return 0
        """

        fact = float(fact)

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._y_scale(fact)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def y_scale_med(self):
        """ @brief Scale y axis by median
        @details Scale the y axis by its median.
        @return 0
        """

        fact = 1/np.median(self.sess.spec.y)

        for s in self.sess.seq:
            try:
                struct = getattr(self.sess, s)
                struct._y_scale(fact)
            except:
                logging.debug(msg_attr_miss(s))
        return 0
