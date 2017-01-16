from astropy import units as u
from astropy.constants import c
from astropy.io import fits
from astropy.table import Column, QTable, Table
import numpy as np
import sys
    
class Spec(QTable):
    """Class for generic spectra"""
    
    ux = u.nm
    uy = u.erg / u.cm ** 2 / u.s / u.nm
    ux_uv_eq = [(u.nm, u.km / u.s, lambda x: np.log(x) *\
        c.to(u.km / u.s).value, lambda x: np.exp(x / c.to(u.km / u.s).value))]
        
    def __init__(self, data = None, masked = None, file = None, bin = None, 
        size = None, flux = None, flux_e = None, ux = None, uy = None):
        """Initialize the spectrum"""
        
        if data is not None:
            QTable.__init__(self, data)
        else:
            try:
                data = fits.open(file)[1].data
                bin = data.field('bin')
                size = data.field('size')
                flux = data.field('flux')            
                flux_e = data.field('flux_e')
            except:
                try:
                    data = fits.open(file)[1].data
                    bin = data.field('WAVEL')
                    size = data.field('PIXSIZE')
                    flux = data.field('FLUX')            
                    flux_e = data.field('FLUXERR')
                except:
                    pass            
                pass            
            
            try:
                bin_unit = bin.unit
            except:
                bin_unit = self.ux
            try:
                size_unit = size.unit
            except:
                size_unit = self.ux
            try:
                flux_unit = flux.unit
            except:
                flux_unit = self.uy
            try:
                flux_e_unit = flux_e.unit
            except:
                flux_e_unit = self.uy
            
            bin_col = Column(bin, name = 'bin', unit = bin_unit)
            size_col = Column(size, name = 'size', unit = size_unit)
            flux_col = Column(flux, name = 'flux', unit = flux_unit)
            flux_e_col = Column(flux_e, name = 'flux_e', unit = flux_e_unit)
            
            QTable.__init__(self, data = (bin_col, size_col, flux_col,
                flux_e_col))
        
    def bin_to_vel(self, cen = 0.0, unit = u.km / u.s, only_size = False):
        """Convert wavelength bins into velocity bins
        
        WRITE TEST!""" 
            
        # Use Table because QTable doesn't accept conversion.
        new = Table(self)
        old_unit = self['bin'].unit
        new['bin'] = new['bin'].to(unit, equivalencies = self.ux_uv_eq)
        new['bin'].unit = unit
        try:
            new['size'].unit = unit
            new['size'][:] = new['bin'][:]
            new['size'][1: len(new) - 1] = 0.5 * (new['size'][2: len(new)] -\
                 new['bin'][0: len(new) - 2])
            new['size'][0] = 1.0 * (new['bin'][1] - new['bin'][0])
            new['size'][len(new) - 1] = 1.0 * (new['bin'][len(new) - 1] -\
                 new['bin'][len(new) - 2])
        except:
            pass
        if only_size is False:
            new['bin'] = new['bin'] - cen.to(unit, equivalencies = 
                self.ux_uv_eq)
        else:
            print(old_unit, self.ux)
            new['bin'] = new['bin'].to(old_unit, equivalencies = 
                self.ux_uv_eq)
            new['bin'].unit = old_unit
        return Spec(new)                
            
    def bin_to_wave(self, cen = None, unit = u.nm, only_size = False):
        """Convert velocity bins into wavelength bins
        
        WRITE TEST!""" 
        
        new = Table(self)
        old_unit = self['bin'].unit
        if only_size is False:
            new['bin'] = new['bin'] + cen.to(old_unit, equivalencies =
                self.ux_uv_eq)
            new['bin'] = new['bin'].to(unit, equivalencies = self.ux_uv_eq)
            new['bin'].unit = unit
        new['size'].unit = unit
        new['size'][:] = new['bin'][:]
        new['size'][1: len(new) - 1] = 0.5 * (new['size'][2: len(new)] -\
             new['bin'][0: len(new) - 2])
        new['size'][0] = 1.0 * (new['bin'][1] - new['bin'][0])
        new['size'][len(new) - 1] = 1.0 * (new['bin'][len(new) - 1] -\
             new['bin'][len(new) - 2])
        return Spec(new)

    def conv_prof(self, prof, mode = 'full'):
        """Convolve a spectrum with a profile
        
        WRITE TEST!"""
    
        if self['bin'][0].si.unit != prof['bin'][0].si.unit:
            print("Units of pixels do not match.")
            sys.exit()
        bin = self['bin'][:] 
        size = self['size'][:]
        flux = np.convolve(self['flux'][:], prof['flux'][:], mode)
        flux_e = np.convolve(self['flux_e'][:], prof['flux'][:], mode)
        shift = int((prof.nrow - 1) * 0.5)
        if len(flux) > len(self): 
            flux = flux[shift: -shift]
            flux_e = flux_e[shift: -shift]
        else: 
            bin = bin[shift: -shift]
            size = size[shift: -shift]
        new = Spec(bin = bin, size = size, flux = flux, flux_e = flux_e,
            ux = self['bin'])
        new.shift = shift
        return new
        
    def save(self, file):
        """
        print(self.columns)
        fits_cols = []
        for col in self.columns:
            print(col, self[col].unit, self[col].value)
            fits_cols.append(fits.Column(name = col, format = 'float32', 
                unit = self[col].unit, array = self[col].value))
        hdu = fits.BinTableHDU.from_columns(fits_cols)
        #hdu.writeto('file')
        """
        new = Table(self)
        new.write(file, overwrite = True)
    
    def to_redshift(self, z = 0):
        """Shift the spectrum to a given redshift"""
        
        self['bin'][:] = self['bin'][:] * (1.0 + z)
        try:
            self['size'][:] = self['size'][:] * (1.0 + z)
        except:
            pass