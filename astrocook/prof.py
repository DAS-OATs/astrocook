from astropy import units as u
from astropy.table import Column, QTable
import astropy.modeling.functional_models as func
import numpy as np
            
class ProfGauss(QTable):
    """Class for Gaussian profiles
    
    A gaussian profile is a a QTable with the following columns:
        -# @bin: channel;
        -# @size: size of the channel; 
        -# @flux: flux density in the channel."""

    def __init__(self, 
                 ampl = 1.0, 
                 stddev = 1e4 * u.m / u.s, 
                 stddev_num = 3, 
                 step = 1e3 * u.m / u.s):
        """Initialize the profile
        
        The profile is defined on a velocity grid with a given @step. The 
        gaussian is centered in zero, with amplitude @ampl and standard
        deviation @stddev. It is cut a @stddev_num from the center.
        Units are initialized from defaults.
        """
        
        stddev_eq = stddev.to(step.unit)
        self.nrow = int(2 * stddev_num * np.ceil(stddev_eq / step) + 1)
        gauss = func.Gaussian1D()
        gauss.amplitude = ampl
        gauss.mean = 0.0
        gauss.stddev = stddev_eq.value
        
        # Data setup
        bin = (np.arange(self.nrow) - (self.nrow - 1) * 0.5) * step.value
        size = np.ones(self.nrow) * step.value
        flux = gauss(bin) / np.sum(gauss(bin))

        # Units setuo
        bin_unit = u.m / u.s
        size_unit = u.m / u.s
        flux_unit = u.erg / u.cm ** 2 / u.s / u.nm
        
        # Column creation
        bin_col = Column(bin, name = 'bin', unit = bin_unit)
        size_col = Column(size, name = 'size', unit = size_unit)
        flux_col = Column(flux, name = 'flux', unit = flux_unit)
        
        # Table creation
        QTable.__init__(self, data = (bin_col, size_col, flux_col))