from .format import Format
from .spectrum import Spectrum
from astropy.io import fits

class Session(object):
    """ Class for sessions.

    A Session is a self-sufficient set of analysis operations."""

    def __init__(self,
                 name):
        self.name = name

    def open(self):
        format = Format()
        hdul = fits.open(self.name)
        try:
            instr = hdul[0].header['INSTRUME']
        except:
            instr = None
            print("INSTRUME not defined")
        try:
            catg = hdul[0].header['HIERARCH ESO PRO CATG']
        except:
            catg = None
            print("HIERARCH ESO PRO CATG not defined")
        try:
            orig = hdul[0].header['ORIGIN']
        except:
            orig = None
            print("ORIGIN not defined")

        # DAS spectrum
        if instr == 'ESPRESSO' and catg == 'RSPEC':
            self.spec = format.das_spectrum(hdul)
        if orig == 'ESO-MIDAS':
            self.spec = format.eso_midas(hdul)
