from .spectrum import Spectrum
from .line_list import LineList
from .syst_list import SystList
from .message import *
from .vars import *
from astropy import units as au
from astropy.table import Table
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
            """
            try:
                meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
            except:
                meta['object'] = ''
                logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
            """
            if struct in ['spec', 'nodes']:
                out = Spectrum(x, xmin, xmax, y, dy, xunit=xunit, yunit=yunit,
                               meta=meta)
            if struct in ['lines']:
                out = LineList(x, xmin, xmax, y, dy, xunit=xunit, yunit=yunit,
                               meta=meta)

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

            if struct in ['lines']:
                try:
                    out._t['source'] = data['source']
                except:
                    pass

            for c in Table(data).colnames:
                if c not in ['x', 'dx', 'xmin', 'xmax', 'y', 'dy']:
                    out._t[c] = data[c]
                    if c in ['conv', 'cont', 'model', 'deabs']:
                        out._t[c].unit = out._t['y'].unit


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
            for c in Table(data).colnames:
                if c not in ['series', 'func', 'z', 'dz', 'logN', 'dlogN', 'b', 'db', 'resol', 'chi2r', 'id', 'z0']:
                    out._t[c] = data[c]

            for k in hdr:
                ks = k.split(' ')
                if 'CONSTR' in ks:
                    if 'ID' in ks: id = int(hdr[k])
                    if 'PAR' in ks: par = hdr[k]
                    if 'VAL' in ks: out._constr['lines_voigt_%i_%s' % (id,par)] \
                        = (id, par, hdr[k])
            #print(out._constr)
        return out

    def eso_adp(self, hdul):
        logging.info(msg_format('ESO Advanced Data Product'))
        """ ESO Advanced Data Product """

        hdr = hdul[0].header
        hdr1 = hdul[1].header
        data = hdul[1].data
        x = data['WAVE'][0]
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['FLUX'][0]
        try:
            dy = data['ERR_FLUX'][0]
        except:
            dy = data['ERR'][0]
        try:
            cont = data['CONTINUUM'][0]
        except:
            cont = []
        xunit = au.Unit(hdr1['TUNIT1']) #au.Angstrom
        yunit = au.Unit(hdr1['TUNIT2']) #au.erg/au.cm**2/au.s/au.Angstrom
        resol = []*len(x)
        meta = hdr #{'instr': hdr['INSTRUME']}
        """
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta, cont=cont)


    def eso_midas_image(self, hdul):
        logging.info(msg_format('ESO MIDAS image'))
        """ ESO-MIDAS image """

        hdr = hdul[0].header
        data = hdul[0].data
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        y = data
        x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]

        xmin, xmax = self._create_xmin_xmax(x)
        dy = np.full(len(x), np.nan)
        resol = []*len(x)
        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        meta = hdr #{'instr': ''}
        """
        try:
            meta['object'] = hdr['FILENAME'][:-4]
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def eso_midas_table(self, hdul):
        logging.info(msg_format('ESO MIDAS table'))
        """ ESO-MIDAS table """

        hdr = hdul[1].header
        data = hdul[1].data
        try:
            x = data['wave']
            xmin = x-data['wpix']*0.5
            xmax = x+data['wpix']*0.5
            y = data['flux']
            dy = data['sigma']
        except:
            x = data['WAVE']
            xmin, xmax = self._create_xmin_xmax(x)
            y = data['NORMFLUX']
            dy = data['STDEV']
        resol = []*len(x)
        #"""
        try:
            cont = data['cont']
        except:
            try:
                cont = data['CONT']
            except:
                if np.all(data['normflux']==y):
                    cont = np.ones(len(x))
                else:
                    cont = []
        #"""
        #cont = []

        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        meta = hdr #{'instr': ''}
        """
        try:
            meta['object'] = hdr['FILENAME'][:-4]
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
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
        meta = hdr #{'instr': 'ESPRESSO'}
        try:
            cont = data['CONT']
        except:
            cont = []
        """
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta, resol=resol,
                        cont=cont)

    def espresso_drs_spectrum(self, hdul):
        """ ESPRESSO DRS S1D format """
        logging.info(msg_format('ESPRESSO DRS S1D'))

        hdr = hdul[0].header
        data = hdul[1].data
        x = data['wavelength']
        xmin, xmax = self._create_xmin_xmax(x)
        try:
            y = data['flux']/(xmax-xmin)#*10#au.nm/au.Angstrom
            dy = data['error']/(xmax-xmin)#*10#au.nm/au.Angstrom
            yunit = au.electron #erg/au.cm**2/au.s/au.nm
        except:
            y = data['flux_cal']/(xmax-xmin)#*10#au.nm/au.Angstrom
            dy = data['error_cal']/(xmax-xmin)#*10#au.nm/au.Angstrom
            yunit = au.erg/au.cm**2/au.s/au.Angstrom
        resol = []*len(x)
        xunit = au.Angstrom
        meta = hdr #{'instr': 'ESPRESSO'}
        """
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
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


    def firehose_spectrum(self, hdul):
        """ FIRE spectrum format, as produced by the FIREhose pipeline  """
        logging.info(msg_format('FIRE'))

        hdr = hdul[0].header
        data = hdul[5].data
        x = data['WAVE'][0]*0.1
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['FLUX'][0]
        dy = data['SIG'][0]
        xunit = au.nm
        yunit = None
        meta = hdr #{'instr': 'FIRE'}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def generic_spectrum(self, hdul):
        """ Generic spectrum """
        logging.info(msg_format('generic'))
        hdr = hdul[0].header
        try:
            zero
        except:
            if len(hdul)>1:
                data = Table(hdul[1].data)
                x_col = np.where([c in data.colnames for c in x_col_names])[0]
                y_col = np.where([c in data.colnames for c in y_col_names])[0]
                dy_col = np.where([c in data.colnames for c in dy_col_names])[0]
                try:
                    x = data[x_col_names[x_col][0]]
                    y = data[y_col_names[y_col][0]]
                    dy = data[dy_col_names[dy_col][0]] if dy_col is not [] \
                        else np.full(len(y), np.nan)
                except:
                    logging.error("I can't recognize columns.")
                    return 0
            else:
                data = hdul[0].data
                x = data[0][:]
                y = data[1][:]
                dy = data[2][:]
            if np.max(x)>3000:
                x = x*0.1
            xmin, xmax = self._create_xmin_xmax(x)
            xunit = au.nm
            yunit = au.erg/au.cm**2/au.s/au.Angstrom
            meta = hdr #{}
            """
            try:
                meta['object'] = hdr['OBJECT']
            except:
                meta['object'] = ''
            """
            return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
        #except:
        #    return None


    def mage_spectrum(self, hdul):
        """ LDSS3 spectrum """
        logging.info(msg_format('QUBRICS'))
        hdr = hdul[0].header
        data = hdul[0].data
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        y = data
        x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        dy = np.full(len(y), np.nan)
        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        meta = hdr #{}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def qubrics_spectrum(self, hdul):
        """ LDSS3 spectrum """
        logging.info(msg_format('QUBRICS'))
        hdr = hdul[0].header
        data = hdul[1].data
        x = data['wave']*0.1
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['flux']
        dy = np.full(len(y), np.nan)
        xunit = au.nm
        yunit = au.erg/au.cm**2/au.s/au.nm
        meta = hdr #{}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def uves_spectrum(self, hdul, hdul_err):
        """ UVES format (two files: FLUXCAL_SCI, FLUXCAL_ERRORBAR_SCI) """
        logging.info(msg_format('uves'))

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        y = hdul[0].data
        dy = hdul_err[0].data
        x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        resol = []*len(x)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom
        meta = hdr #{'instr': 'UVES'}
        meta['v_bary'] = hdr['HIERARCH ESO QC VRAD BARYCOR']
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('OBJECT'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def uves_popler_spectrum(self, hdul):
        """ UVES_popler format """
        logging.info(msg_format('UVES_popler'))

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        data = hdul[0].data
        y = data[:][0]*data[:][3]
        #dy = data[:][1]#*data[:][3]
        dy = data[:][2]*data[:][3]
        x = 10**np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        resol = []*len(x)
        cont = data[:][3]#np.ones(len(x))
        sel = np.where(y != 1)
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom
        meta = hdr #{'instr': 'UVES'}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('OBJECT'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

    def wfccd_spectrum(self, hdul):
        """ WFCCD format """
        logging.info(msg_format('WFCCD'))

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CD1_1']
        naxis1 = hdr['NAXIS1']
        y = hdul[0].data[0][0]
        dy = hdul[0].data[1][0]
        x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.nm
        meta = hdr #{'instr': 'WFCCD'}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('OBJECT'))
        """
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
        xunit = au.Angstrom
        yunit = au.electron/au.Angstrom #erg/au.cm**2/au.s/au.nm
        meta = hdr #{'instr': 'X-shooter'}
        """
        try:
            meta['object'] = hdul._file.name.split('/')[-1][:-5]
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
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
        meta = hdr #{'instr': 'X-shooter'}
        """
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def xshooter_vacbary_spectrum(self, hdul):
        """ X-shooter VAC/BARY format """
        logging.info(msg_format('X-shooter VAC/BARY'))

        logging.info(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        hdr = hdul[0].header
        data = hdul[1].data
        x = data['wave']
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['flux']
        dy = data['err']
        resol = []*len(x)
        xunit = au.nm
        yunit = au.electron/au.nm #erg/au.cm**2/au.s/au.nm
        meta = hdr #{'instr': 'X-shooter'}
        """
        try:
            meta['object'] = hdr['HIERARCH ESO OBS TARG NAME']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)


    def xshooter_merge1d_spectrum(self, hdul):
        logging.info(msg_format('X-shooter MERGE1D'))

        hdr = hdul[0].header
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        y = hdul[0].data
        dy = hdul[1].data
        x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)]
        xmin, xmax = self._create_xmin_xmax(x)
        resol = []*len(x)
        xunit = au.nm
        yunit = au.erg/au.cm**2/au.s/au.nm
        meta = hdr #{'instr': 'X-shooter'}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('HIERARCH ESO OBS TARG NAME'))
        """
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
        meta = hdr #{'instr': 'X-shooter'}
        """
        try:
            meta['object'] = hdr['OBJECT']
        except:
            meta['object'] = ''
            logging.warning(msg_descr_miss('OBJECT'))
        """
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
