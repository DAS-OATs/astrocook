from .cookbook_flux_old import CookbookFluxOld

from astropy import units as au
import numpy as np

class CookbookFlux(CookbookFluxOld):
    """ Cookbook of utilities for flux calibration
    @details This cookbook contains utilities to rescale the flux, correct
    it for reddening and Lyman-alpha opacity, and computing the flux
    cross-correlation function.
    """

    def __init__(self):
        super(CookbookFlux, self).__init__()


    def rebin(self, xstart=None, xend=None, dx=10.0, xunit=au.km/au.s,
              kappa=None, norm=False, filling=np.nan):
        """ @brief Re-bin spectrum
        @details Apply a new binning to a spectrum, with a constant bin size.
        @url flux_cb.html#re-bin-spectrum
        @param xstart Start wavelength (nm)
        @param xend End wavelength (nm)
        @param dx Step in x
        @param xunit Unit of wavelength or velocity
        @param kappa Number of sigma to clip outliers
        @param norm Return normalized spectrum, if continuum exists
        @param filling Value to fill region without data
        @return Session with rebinned spectrum
        """

        try:
            xstart = None if xstart in [None, 'None'] else float(xstart)
            xend = None if xend in [None, 'None'] else float(xend)
            dx = float(dx)
            xunit = au.Unit(xunit)
            kappa = None if kappa in [None, 'None'] else float(kappa)
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

        spec_out = spec_in._rebin(xstart, xend, dx, xunit, y, dy, kappa, filling)
        if cont:
            if not norm:
                try:  # x-axis in wavelengths
                    spec_out.t['cont'] = np.interp(
                        spec_out.x.to(au.nm).value, spec_in.x.to(au.nm).value,
                        spec_in.t['cont']) * spec_out.y.unit
                except:  # x-axis in velocities
                    spec_out.t['cont'] = np.interp(
                        spec_out.x.to(au.km/au.s).value,
                        spec_in.x.to(au.km/au.s).value,
                        spec_in.t['cont']) * spec_out.y.unit
            else:
                spec_out.t['cont'] = np.ones(len(spec_out.t)) * spec_out.y.unit

        # Create a new session
        from .session import Session
        new = Session(gui=self.sess._gui, name=self.sess.name+'_rebinned',
                      spec=spec_out)
        return new


    def smooth(self):
        """@brief Smooth spectrum 🚧
        @details 🚧
        @url flux_cb.html#smooth-spectrum
        """

        return 0


    def rescale(self):
        """@brief Rescale spectrum 🚧
        @details 🚧
        @url flux_cb.html#rescale-spectrum
        """

        return 0


    def deredden(self, ebv=0.03, rv=3.1):
        """@brief De-redden spectrum
        @details Correct the spectrum flux for reddening due to extinction.
        @url flux_cb.html#de-redden-spectrum
        @param ebv Color excess
        @param rv Ratio of total selective extinction
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


    def adjust_mags(self, bands, refs, deg=1):
        """ @brief Adjust magnitudes
        @details Adjust the flux of the spectrum to its AB magnitudes.
        @url flux_cb.html#adjust-magnitudes
        @param bands List of bands ('u', 'v', 'g', 'r', 'i', or 'z')
        @param refs List of reference magnitudes
        @param deg Degree of polynomial regression
        @return 0
        """

        mags_in = self._mags_compute(bands)
        x = np.array([self._filt_x[b] for b in bands])
        diffs = np.array([float(refs[b])-mags_in[b] for b in bands])

        x = x[~np.isnan(diffs)]
        diffs = diffs[~np.isnan(diffs)]

        spec = self.sess.spec
        corr = 10**(np.poly1d(np.polyfit(x, diffs, deg)/-2.5)\
                (spec.x.to(au.nm).value))
        spec._t['y'] = spec._t['y']*corr
        spec._t['dy'] = spec._t['dy']*corr
        mags_in = self._mags_compute(bands)

        return 0


    def estimate_snr(self):
        """ @brief Estimate SNR
        @details Estimate the signal-to-noise ratio per pixel.
        @url flux_cb.html#estimate-snr
        @return 0
        """

        spec = self.sess.spec
        if 'snr' not in spec._t.colnames:
            logging.info("I'm adding column 'snr'.")
        else:
            logging.info("I'm updating column 'snr'.")

        spec._t['snr'] = spec.y/spec.dy

        return 0


    def estimate_rms(self, hwindow=10000, std=20.0):
        """ @brief Estimate RMS
        @details Estimate flux error by computing the root-mean-square (RMS) of
        the flux within a running window.
        @url flux_cb.html#estimate-rms
        @param hwindow Half-size in pixels of the running window
        @param std Standard deviation of the gaussian (km/s)
        @return 0
        """

        try:
            hwindow = int(hwindow)
            std = float(std) * au.km/au.s
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
        spec._t['y_rms'] = at.Column(y_rms, dtype=float)*spec._t['dy'].unit
        self.sess.spec._gauss_convolve(std, 'y_rms', 'dy')

        return 0
