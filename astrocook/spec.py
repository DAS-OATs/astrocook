from astropy import units as u
from astropy.constants import c
from astropy.io import fits
from astropy.table import Column, QTable
import copy
import numpy as np
import sys
    
class Spec(QTable):
    """Class for generic spectra
    
    A generic spectrum is a QTable with the following columns: 
        -# @exp: exposure;
        -# @ord: spectral order;
        -# @bin: channel;
        -# @size: size of the channel; 
        -# @qual: quality of the channel;
        -# @resol: spectral resolution in the channel; 
        -# @resol_e: error on @resol; 
        -# @flux: flux density in the channel;
        -# @flux_e: error on @flux.
        
    The following metadata are also expected:
        -# @nexp: number of exposures;
        -# @exptime: integration time of all exposures;
        -# @misc: miscellaneous information (TBD).        
    """
    
    def __init__(self, 
                 data = None, 
                 masked = None, 
                 file = None, 
                 exptime = None,
                 misc = None):
        """Initialize the spectrum
        
        The spectrum can be initialized either from @data (default) or from 
        @file. In the latter case, a FITS spectrum is expected in input, either
        from astrocook or from the ESPRESSO DAS. 
        Units are initialized either from @data/@file or from defaults.
        Metadata are defined as properties (see methods below).
        
        TBD: 
            1. A method should be written to cope with different formats.
            2. The units setup could be rewritten in a more clever way.
            3. The number of spectral orders for all exposures should be 
               computed.
        """
        
        if data is not None:
            QTable.__init__(self, data)
        else:
            
            # From astrocook
            try:
                data = fits.open(file)[1].data
                exp = data.field('exp')
                ord = data.field('ord')
                bin = data.field('bin')
                size = data.field('size')
                qual = data.field('qual')
                resol = data.field('resol')
                resol_e = data.field('resol_e')
                flux = data.field('flux')            
                flux_e = data.field('flux_e')
            except:
                
                # From the ESPRESSO DAS
                try:
                    data = fits.open(file)[1].data
                    exp = [1] * len(data)
                    ord = [0] * len(data)
                    bin = data.field('WAVEL')
                    size = data.field('PIXSIZE')
                    qual = [1] * len(data)
                    resol = [60000.] * len(data)
                    resol_e = [1000.] * len(data)
                    flux = data.field('FLUX')            
                    flux_e = data.field('FLUXERR')
                except:
                    sys.exit("FITS file not recognized.")
            
            # Units setup
            try:
                bin_unit = bin.unit
            except:
                bin_unit = u.nm
            try:
                size_unit = size.unit
            except:
                size_unit = u.nm
            try:
                flux_unit = flux.unit
            except:
                flux_unit = u.erg / u.cm ** 2 / u.s / u.nm
            try:
                flux_e_unit = flux_e.unit
            except:
                flux_e_unit = u.erg / u.cm ** 2 / u.s / u.nm
            
            # Column creation
            exp_col = Column(exp, name = 'exp')
            ord_col = Column(ord, name = 'ord')
            bin_col = Column(bin, name = 'bin', unit = bin_unit)
            size_col = Column(size, name = 'size', unit = size_unit)
            qual_col = Column(qual, name = 'qual')
            resol_col = Column(resol, name = 'resol')
            resol_e_col = Column(resol_e, name = 'resol_e')
            flux_col = Column(flux, name = 'flux', unit = flux_unit)
            flux_e_col = Column(flux_e, name = 'flux_e', unit = flux_e_unit)
            
            # Table creation
            QTable.__init__(self, data = (exp_col, ord_col, bin_col, size_col, 
                    qual_col, resol_col, resol_e_col, flux_col, flux_e_col))

        # Metadata
        self._nexp = np.max(self['exp'])
        if exptime == None:
            exptime = []
        elif (len(exptime) != self._nexp):
            print("Exposure time array has a wrong size and will be "
                  "neglected.")
            exptime = []
        self._exptime = copy.deepcopy(exptime)        
        if misc == None:
            misc = {}
        self._misc = copy.deepcopy(misc)

    @property
    def nexp(self):
        """Property: number of exposures"""

        return self._nexp

    @property
    def exptime(self):
        """Property: integration time of all exposures"""

        return self._exptime
        
    @exptime.setter
    def exptime(self, value):
        """@exptime setter"""

        if value and len(value) == self._nexp:
            self._exptime = value
        else:
            print("Wrong number of exposure times: neglected.") 

    @property
    def misc(self):
        """Property: miscellaneous information"""
        
        return self._misc
        
    @misc.setter
    def misc(self, value):
        """@misc setter"""

        self._misc = self._misc.update({value})