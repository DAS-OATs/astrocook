from astropy import units as au
from .spectrum import Spectrum

class Format(object):
    """ Class for file formats. """

    def __init__(self):
        pass

    def das_spectrum(self, hdul):
        """ DAS FSPEC/RSPEC format """

        hdr = hdul[0].header
        data = hdul[1].data
        x = data['WAVEL']
        xmin = data['WAVEL']-data['PIXSIZE']*0.5
        xmax = data['WAVEL']+data['PIXSIZE']*0.5
        y = data['FLUX']
        dy = data['FLUXERR']
        xunit = au.nm
        yunit = au.electron #erg/au.cm**2/au.s/au.nm
        try:
            meta = {'object': hdr['HIERARCH ESO OBS TARG NAME']}
        except:
            meta = {'object': ""}
            print("HIERARCH ESO OBS TARG NAME not defined")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def eso_midas(self, hdul):
        """ ESO-MIDAS Table """

        hdr = hdul[1].header
        data = hdul[1].data
        x = data[hdr['TTYPE1']]
        xmin = data[hdr['TTYPE1']]
        xmax = data[hdr['TTYPE1']]
        y = data[hdr['TTYPE2']]
        dy = data[hdr['TTYPE3']]
        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        try:
            meta = {'object': hdr['HIERARCH ESO OBS TARG NAME']}
        except:
            meta = {'object': ""}
            print("HIERARCH ESO OBS TARG NAME not defined")
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
