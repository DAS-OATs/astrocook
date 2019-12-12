from . import version
from .cookbook import Cookbook
from .format import Format
from .functions import detect_local_minima
from .line_list import LineList
from .message import *
#from .model import Model
from .spectrum import Spectrum
from .syst_list import SystList
from .syst_model import SystModel
#from .model_list import ModelList
from .vars import *
#from astropy import constants as ac
from astropy import units as au
from astropy.io import ascii, fits
from copy import deepcopy as dc
import logging
from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.signal import argrelmin
import tarfile
import time

class Session(object):
    """ Class for sessions.

    A Session is a self-sufficient set of analysis operations."""

    def __init__(self,
                 path=None,
                 name=None,
                 spec=None,
                 spec_form=None,
                 nodes=None,
                 lines=None,
                 systs=None,
                 mods=None,
                 twin=False):
        #self._gui = gui
        self.path = path
        self.name = name
        self.spec = spec
        self.spec_form = spec_form
        self.nodes = nodes
        self.lines = lines
        self.systs = systs
        self.mods = mods
        self.seq = seq  # From .vars
        self.cb = Cookbook(self)
        self._open_twin = twin


    def _append(self, frame, append=True):
        if append and hasattr(self, frame.__name__):
            getattr(self, frame.__name__)._append(frame)
        else:
            setattr(self, frame.__name__, frame)

    def open(self):

        format = Format()
        if self.path[-3:] == 'acs':
            root = '/'.join(self.path.split('/')[:-1])
            with tarfile.open(self.path) as arch:
                arch.extractall(path=root)
                hdul = fits.open(self.path[:-4]+'_spec.fits')
                hdr = hdul[1].header
        else:
            hdul = fits.open(self.path)
            hdr = hdul[0].header

        try:
            instr = hdr['INSTRUME']
        except:
            instr = 'undefined'
        try:
            orig = hdr['ORIGIN']
        except:
            orig = 'undefined'
        try:
            catg = hdr['HIERARCH ESO PRO CATG']
        except:
            catg = 'undefined'

        try:
            hist = [i.split(' ') for i in str(hdr['HISTORY']).split('\n')]
            hist = [i for j in hist for i in j]
            if 'UVES_popler:' in hist:
                instr = 'UVES'
                orig = 'POPLER'
            if 'XSHOOTER_REDUCE' in hist:
                instr = 'XSHOOTER'
                orig = 'REDUCE'
        except:
            pass

        try:  # For XQR-30 spectra
            ttype = [hdul[1].header['TTYPE%i' % i].strip() for i in range(1,6)]
            if ttype[1] == 'FLUX' and ttype[3] == 'FLUX_NOCORR':
                instr = 'XSHOOTER'
                orig = 'XQR-30'
        except:
            pass

        if instr == None:
            logging.warning(msg_descr_miss('INSTRUME'))
        if catg == None:
            logging.warning(msg_descr_miss('HIERARCH ESO PRO CATG'))
        if orig == None:
            logging.warning(msg_descr_miss('ORIGIN'))

        # Astrocook structures
        logging.debug("Instrument: %s; origin: %s; category: %s."
                      % (instr, orig, catg))
        if orig[:9] == 'Astrocook':
            for s in self.seq:
                try:
                    hdul = fits.open(self.path[:-4]+'_'+s+'.fits')
                    setattr(self, s, format.astrocook(hdul, s))
                    os.remove(self.path[:-4]+'_'+s+'.fits')
                except:
                    pass
            if self.spec is not None and self.systs is not None:
                self.cb._mods_recreate()

        # ESO-MIDAS spectrum
        if orig == 'ESO-MIDAS':
            self.spec = format.eso_midas(hdul)

        # ESPRESSO DRS spectrum
        if instr == 'ESPRESSO' and catg[0:3] == 'S1D':
            self.spec = format.espresso_drs_spectrum(hdul)
            self.spec_form = format.espresso_spectrum_format(
                ascii.read('espr_spec_form.dat'))

        # ESPRESSO DAS spectrum
        if instr == 'ESPRESSO' and catg[1:5] == 'SPEC':
            self.spec = format.espresso_das_spectrum(hdul)
            self.spec_form = format.espresso_spectrum_format(
                ascii.read('espr_spec_form.dat'))

        # UVES POPLER spectrum
        if instr == 'UVES' and orig == 'POPLER':
            self.spec = format.uves_popler_spectrum(hdul)

        # XQR-30 spectrum
        if instr == 'XSHOOTER' and orig == 'XQR-30':
            self.spec = format.xqr30_spectrum(hdul, corr=self._open_twin)
            self._open_twin = not self._open_twin

        # XSHOOTER DAS spectrum
        if instr == 'XSHOOTER' and catg[1:5] == 'SPEC':
            self.spec = format.xshooter_das_spectrum(hdul)

        # XSHOOTER_REDUCE spectrum
        if instr == 'XSHOOTER' and orig == 'REDUCE':
            hdul_e = fits.open(self.path[:-5]+'e.fits')
            self.spec = format.xshooter_reduce_spectrum(hdul, hdul_e)


    def save(self, path):

        root = path[:-4]
        stem = root.split('/')[-1]
        with tarfile.open(root+'.acs', 'w:gz') as arch:
            for s in self.seq:
                if hasattr(self, s) and getattr(self, s) is not None:
                    if s=='systs':
                        try:
                            np.savetxt(root+'_compl.dat', self.compl, fmt='%s')
                        except:
                            pass
                        try:
                            np.savetxt(root+'_corr.dat', self.corr, fmt='%s')
                        except:
                            pass
                        try:
                            np.savetxt(root+'_merge.dat', self.merge, fmt='%s')
                        except:
                            pass
                    name = root+'_'+s+'.fits'
                    obj = dc(getattr(self, s))
                    t = obj._t
                    if s!='systs':
                        t['x'] = t['x'].to(au.nm)
                        t['xmin'] = t['xmin'].to(au.nm)
                        t['xmax'] = t['xmax'].to(au.nm)
                    t.meta = obj._meta
                    t.meta['ORIGIN'] = 'Astrocook'
                    #t.meta['HIERARCH ASTROCOOK VERSION'] = version
                    #t.meta['HIERARCH ASTROCOOK STRUCT'] = s
                    for c in t.colnames:
                        t[c].unit = au.dimensionless_unscaled
                    t.write(name, format='fits', overwrite=True)
                    arch.add(name, arcname=stem+'_'+s+'.fits')
                    os.remove(name)
                    logging.info("I've saved frame %s as %s."
                                 % (s, stem+'_'+s+'.fits'))
                else:
                    logging.error("I haven't found any frame %s to save." % s)
