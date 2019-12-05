from .spectrum import Spectrum
from .line_list import LineList
from .syst_list import SystList
from .message import *
from astropy import units as au
import logging
import numpy as np

class Format(object):
    """ Class for file formats. """

    def __init__(self):
        pass

    def _create_xmin_xmax(self, x):
        mean = 0.5*(x[1:]+x[:-1])
        xmin = np.append(x[0], mean)
        xmax = np.append(mean, x[-1])
        return xmin, xmax

    def astrocook(self, hdul, struct):
        logging.info(msg_format('Astrocook %s' % struct))
        hdr = hdul[1].header
        data = hdul[1].data

        if struct in ['spec', 'nodes', 'lines']:

            x = data['x']
            xmin = data['xmin']
            xmax = data['xmax']
            y = data['y']
            dy = data['dy']
            xunit = au.nm
            yunit = au.erg/au.cm**2/au.s/au.Angstrom
            meta = hdr
            try:
                meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
            except:
                meta['object'] = ''
                logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
            if struct in ['spec', 'nodes']:
                out = Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
            if struct in ['lines']:
                out = LineList(x, xmin, xmax, y, dy, xunit, yunit, meta)

            # Additional columns
            if struct in ['spec']:
                try:
                    out._t['cont'] = data['cont']
                    out._t['cont'].unit = out._t['y'].unit
                except:
                    pass
                try:
                    out._t['conv'] = data['conv']
                    out._t['conv'].unit = out._t['y'].unit
                except:
                    pass
                try:
                    out._t['line_mask'] = data['line_mask']
                except:
                    pass
                try:
                    out._t['model'] = data['model']
                    out._t['model'].unit = out._t['y'].unit
                except:
                    pass
                try:
                    out._t['deabs'] = data['deabs']
                    out._t['deabs'].unit = out._t['y'].unit
                except:
                    pass
                try:
                    out._t['resol'] = data['resol']
                except:
                    pass

        if struct in ['systs']:
            series = data['series']
            func = data['func']
            z = data['z']
            dz = data['dz']
            logN = data['logN']
            dlogN = data['dlogN']
            b = data['b']
            db = data['db']
            resol = data['resol']
            chi2r = data['chi2r']
            id = data['id']
            out = SystList(func=func, series=series, z=z, dz=dz, logN=logN,
                           dlogN=dlogN, b=b, db=db, resol=resol, chi2r=chi2r,
                           id=id)
            out._t['z0'] = data['z0']

        return out

    def eso_midas(self, hdul):
        logging.info(msg_format('ESO MIDAS'))
        """ ESO-MIDAS Table """

        hdr = hdul[1].header
        data = hdul[1].data
        x = data['wave']
        xmin = x-data['wpix']*0.5
        xmax = x+data['wpix']*0.5
        y = data['flux']
        dy = data['sigma']
        resol = []*len(x)
        #"""
        try:
            cont = data['cont']
        except:
            try:
                cont = data['CONT']
            except:
                cont = []
        #"""
        #cont = []

        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        meta = {'instr': ''}
        try:
            meta['object'] = hdr['FILENAME'][:-4]
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta, cont=cont)

    def espresso_das_spectrum(self, hdul):
        """ ESPRESSO DAS FSPEC/RSPEC format """
        logging.info(msg_format('ESPRESSO DAS'))

        hdr = hdul[0].header
        data = hdul[1].data
        x = data['WAVEL']
        xmin = x-data['PIXSIZE']*0.5
        xmax = x+data['PIXSIZE']*0.5
        y = data['FLUX']
        dy = data['FLUXERR']
        resol = np.array([70000]*len(x))
        xunit = au.nm
        yunit = au.electron/au.nm #erg/au.cm**2/au.s/au.nm
        meta = {'instr': 'ESPRESSO'}
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta, resol=resol)

    def espresso_drs_spectrum(self, hdul):
        """ ESPRESSO DRS S1D format """
        logging.info(msg_format('ESPRESSO DRS S1D'))

        hdr = hdul[0].header
        data = hdul[1].data
        x = data['wavelength']
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['flux']/(xmax-xmin)#*10#au.nm/au.Angstrom
        dy = data['error']/(xmax-xmin)#*10#au.nm/au.Angstrom
        resol = []*len(x)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom #erg/au.cm**2/au.s/au.nm
        meta = {'instr': 'ESPRESSO'}
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def espresso_spectrum_format(self, data):
        """ ESPRESSO spectrum format """
        logging.info(msg_format('ESPRESSO'))

        x = data['col2']
        xmin = data['col8']
        xmax = data['col9']
        y = np.zeros(len(x))
        dy = np.zeros(len(x))
        resol = []*len(x)
        xunit = au.nm
        yunit = None
        meta = {}
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def uves_popler_spectrum(self, hdul):
        """ UVES_popler format """
        logging.info(msg_format('UVES_popler'))

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        data = hdul[0].data
        y = data[:][0]#*data[:][3]
        #dy = data[:][1]#*data[:][3]
        dy = data[:][2]#*data[:][3]
        x = 10**np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        resol = []*len(x)
        cont = np.ones(len(x))
        sel = np.where(y != 1)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom
        meta = {'instr': 'UVES'}
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('OBJECT'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def xshooter_das_spectrum(self, hdul):
        """ X-shooter DAS FSPEC/RSPEC format """
        logging.info(msg_format('X-shooter DAS'))

        logging.info(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        hdr = hdul[0].header
        data = hdul[1].data
        x = data['WAVEL']
        xmin = x-data['PIXSIZE']*0.5
        xmax = x+data['PIXSIZE']*0.5
        y = data['FLUX']
        dy = data['FLUXERR']
        resol = []*len(x)
        xunit = au.nm
        yunit = au.electron/au.nm #erg/au.cm**2/au.s/au.nm
        meta = {'instr': 'X-shooter'}
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def xqr30_spectrum(self, hdul, corr=True):
        logging.info(msg_format('xqr30'))

        hdr = hdul[1].header
        data = hdul[1].data
        x = data['WAVE'][0]
        xmin, xmax = self._create_xmin_xmax(x)
        if corr:
            y = data['FLUX'][0]
            dy = data['ERROR'][0]
        else:
            y = data['FLUX_NOCORR'][0]
            dy = data['ERROR_NOCORR'][0]
        resol = []*len(x)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom #erg/au.cm**2/au.s/au.nm
        meta = {'instr': 'X-shooter'}
        try:
            meta['object'] = hdul._file.name.split('/')[-1][:-5]
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def xshooter_reduce_spectrum(self, hdul, hdul_e):
        logging.info(msg_format('xshooter_reduce'))

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        data = hdul[0].data
        data_e = hdul_e[0].data
        y = data
        dy = data_e
        x = 10**np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        resol = []*len(x)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom
        meta = {'instr': 'X-shooter'}
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('OBJECT'))
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
