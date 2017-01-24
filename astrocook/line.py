from astropy import units as u
from astropy.io import fits
from astropy.table import Column, QTable
import copy
import sys

class Line(QTable):
    """Class for generic line lists
    
    A generic line list is a QTable with the following columns: 
        -# @bin: channel;
        -# @qual: quality of the channel;
        -# @flux: flux density in the channel;
        -# @flux_e: error on @flux.

    The following metadata are also expected:
        -# @misc: miscellaneous information (TBD).        
    """

    def __init__(self, 
                 data = None, 
                 masked = None, 
                 file = None, 
                 misc = None):
        """Initialize the line list
        
        The spectrum can be initialized either from @data (default) or from 
        @file. In the latter case, a FITS spectrum is expected in input, either
        from astrocook or from the ESPRESSO DAS. 
        Units are initialized either from @data/@file or from defaults.
        Metadata are defined as properties (see methods below).
        
        TBD: 
            1. A method should be written to cope with different formats.
            2. The units setup could be rewritten in a more clever way.
        """

        if data is not None:
            QTable.__init__(self, data)
        else:
        
            # From astrocook
            try:
                data = fits.open(file)[1].data
                bin = data.field('bin')
                qual = data.field('qual')
                flux = data.field('flux')            
                flux_e = data.field('flux_e')
            except:
            
                # From the ESPRESSO DAS
                try:
                    data = fits.open(file)[1].data
                    bin = data.field('WAVEL')
                    qual = [1] * len(data)
                    flux = data.field('PEAK')            
                    flux_e = data.field('PEAK')
                except:
                    sys.exit("FITS file not recognized.")
            
            # Units setup
            try:
                bin_unit = bin.unit
            except:
                bin_unit = u.nm
            try:
                flux_unit = flux.unit
            except:
                flux_unit = u.erg / u.cm ** 2 / u.s / u.nm
            try:
                flux_e_unit = flux_e.unit
            except:
                flux_e_unit = u.erg / u.cm ** 2 / u.s / u.nm
            
            # Column creation
            bin_col = Column(bin, name = 'bin', unit = bin_unit)
            qual_col = Column(qual, name = 'qual')
            flux_col = Column(flux, name = 'flux', unit = flux_unit)
            flux_e_col = Column(flux_e, name = 'flux_e', unit = flux_e_unit)

            # Table creation
            QTable.__init__(self, data = (bin_col, qual_col, flux_col, 
                    flux_e_col))
        
        # Metadata
        if misc == None:
            misc = {}
        self._misc = copy.deepcopy(misc)

    @property
    def misc(self):
        """Property: miscellaneous information"""
        
        return self._misc
        
    @misc.setter
    def misc(self, value):
        """@misc setter"""

        self._misc = self._misc.update({value})