from .vars import *
from .message import *
from copy import deepcopy as dc
import numpy as np
import sys
from tqdm import tqdm

class CookbookGeneral(object):
    """ Cookbook of general utilities
    """

    def __init__(self):
        pass

    def gauss_convolve(self, std=5, input_col='y', output_col='conv'):
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


    def region_extract(self, xmin, xmax):
        """ @brief Extract region
        @details Extract a spectral region as a new frame.
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
                logging.debug("Attribute %s does not support region"
                              "extraction." % s)
                try:
                    kwargs[s] = getattr(self.sess, s)
                except:
                    logging.debug(msg_attr_miss(s))
                    kwargs[s] = None
        if kwargs['spec'] != None:
            from .session import Session
            new = Session(**kwargs)
        else:
            new = None

        return new


    def rebin(self, dx=10.0, xunit=au.km/au.s):
        """ @brief Rebin spectrum
        @details Rebin a spectrum with a given velocity step. A new session is
        created with the rebinned spectrum. Other objects from the old session
        (line lists, etc.) are discarded.
        @param dx Step in x
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        try:
            dx = float(dx)
            xunit = au.Unit(xunit)
        except ValueError:
            logging.error(msg_param_fail)
            return None

        # A deep copy is created, so the original spectrum is preserved
        spec_in = dc(self.sess.spec)
        spec_out = spec_in._rebin(dx, xunit)

        # Create a new session
        from .session import Session
        new = Session(name=self.sess.name+'_rebinned', spec=spec_out)
        return new


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


    def x_convert(self, zem=0, xunit=au.km/au.s):
        """ @brief Convert x axis
        @details Convert the x axis to wavelength or velocity units.
        @param zem Emission redshift, to use as a 0-point for velocities
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
        @details Convert the y axis to flux density units.
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
