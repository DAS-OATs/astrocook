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

    def astrocook(self, frame, struct):
        logging.info(msg_format('Astrocook %s' % struct))
        try:
            hdr = frame[1].header
            data = frame[1].data
        except:
            hdr = {}
            data = frame


        if struct in ['spec', 'nodes', 'lines']:

            x = data['x']
            xmin = data['xmin']
            xmax = data['xmax']
            y = data['y']
            dy = data['dy']
            xunit = au.nm
            yunit = au.erg/au.cm**2/au.s/au.Angstrom
            meta = hdr[:50]
            delete = []
            for m in meta:
                if 'AC' not in m and len(m)>8:
                    delete.append(m)
            for d in delete:
                del meta[d]
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
            z0 = data['z0']
            z = data['z']
            dz = data['dz']
            logN = data['logN']
            dlogN = data['dlogN']
            b = data['b']
            db = data['db']
            try:
                btur = data['btur']
                dbtur = data['dbtur']
            except:
                btur = np.zeros(len(data))
                dbtur = [np.nan]*len(data)
            resol = data['resol']
            chi2r = data['chi2r']
            id = data['id']
            out = SystList(func=func, series=series, z0=z0, z=z, dz=dz,
                           logN=logN, dlogN=dlogN, b=b, db=db, btur=btur,
                           dbtur=dbtur, resol=resol, chi2r=chi2r, id=id)
            out._t.sort('z')
            for c in Table(data).colnames:
                if c not in ['series', 'func', 'z', 'dz', 'logN', 'dlogN', 'b',
                             'db', 'resol', 'chi2r', 'id', 'z0']:
                    out._t[c] = data[c]
            for k in hdr:
                ks = k.split(' ')
                if 'CONSTR' in ks:
                    id_check = 'AC CONSTR ID '+ks[-1]
                    if 'id' in out._t.colnames and hdr[id_check] in out._t['id']:
                        if 'ID' in ks: id = int(hdr[k])
                        if 'PAR' in ks: par = hdr[k]
                        if 'VAL' in ks: out._constr['lines_voigt_%i_%s' % (id,par)] \
                            = (id, par, hdr[k])
        return out


    def eso_adp(self, hdul):
        logging.info(msg_format('ESO Advanced Data Product'))
        """ ESO Advanced Data Product """

        hdr = hdul[0].header
        hdr1 = hdul[1].header
        data = Table(hdul[1].data)
        x = data['WAVE'][0]
        xmin, xmax = self._create_xmin_xmax(x)
        y = data['FLUX'][0]
        try:
            dy = data['ERR_FLUX'][0]
        except:
            dy = data['ERR'][0]

        #"""
        try:
            cont = data['CONTINUUM'][0]
        except:
            cont = []
        #"""
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
        spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta, cont=cont)
        #spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

        for i,c in enumerate(data.colnames):
            if c not in ['WAVE', 'FLUX', 'ERR_FLUX', 'ERR', 'CONTINUUM']:
                spec._t[c] = data[c][0]
                spec._t[c].unit = hdr1['TUNIT%i' % (i+1)]

        return spec

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
        dy = 0.05*y #np.full(len(x), np.nan)
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
        data = Table(hdul[1].data)

        x_col_names = np.array(['wave', 'WAVE'])
        y_col_names = np.array(['fluxc', 'flux', 'FLUX'])
        dy_col_names = np.array(['errc', 'err', 'sigma', 'STDEV'])
        wpix_col_names = np.array(['wpix', 'WPIX'])
        cont_col_names = np.array(['cont', 'CONT', 'Continuum'])
        norm_col_names = np.array(['normflux', 'NORMFLUX'])

        x_col = np.where([c in data.colnames for c in x_col_names])[0]
        y_col = np.where([c in data.colnames for c in y_col_names])[0]
        dy_col = np.where([c in data.colnames for c in dy_col_names])[0]
        wpix_col = np.where([c in data.colnames for c in wpix_col_names])[0]
        cont_col = np.where([c in data.colnames for c in cont_col_names])[0]
        norm_col = np.where([c in data.colnames for c in norm_col_names])[0]

        try:
            x = data[x_col_names[x_col][0]]
            y = data[y_col_names[y_col][0]]
            dy = data[dy_col_names[dy_col][0]] if dy_col!=[] \
                else np.full(len(y), np.nan)
        except:
            logging.error("I can't recognize columns.")
            return 0

        try:
            xmin = x-data[wpix_col_names[wpix_col][0]]*0.5
            xmax = x+data[wpix_col_names[wpix_col][0]]*0.5
        except:
            xmin, xmax = self._create_xmin_xmax(x)

        try:
            cont = data[cont_col_names[cont_col][0]]
        except:
            try:
                if np.all(data[norm_col_names[norm_col][0]]==y):
                    cont = np.ones(len(x))
            except:
                cont = []
        """
        try:
            x = data['wave']
            try:
                xmin = x-data['wpix']*0.5
                xmax = x+data['wpix']*0.5
            except:
                xmin = x-data['Wpix']*0.5
                xmax = x+data['Wpix']*0.5
            y = data['flux']
            try:
                dy = data['sigma']
            except:
                dy = data['err']
        except:
            x = data['WAVE']
            xmin, xmax = self._create_xmin_xmax(x)
            y = data['NORMFLUX']
            dy = data['STDEV']
        resol = []*len(x)

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

        """

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

    def espresso_s1d_spectrum(self, hdul):
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


    def espresso_s2d_spectrum(self, hdul, row=None, slice=None):


        """ ESPRESSO DRS S2D format """
        logging.info(msg_format('ESPRESSO DRS S2D'))
        hdr = hdul[0].header

        span = 0
        if row is None and slice is None:
            x = np.ravel(hdul['WAVEDATA_VAC_BARY'].data[:,span:-span-1])
            y = np.ravel(hdul['SCIDATA'].data[:,span:-span-1])
            dy = np.ravel(hdul['ERRDATA'].data[:,span:-span-1])
            q = np.ravel(hdul['QUALDATA'].data[:,span:-span-1])
        elif row is None:
            r = range(slice, hdul['WAVEDATA_VAC_BARY'].data.shape[0], 2)
            x = np.ravel(hdul['WAVEDATA_VAC_BARY'].data[r,span:-span-1])
            y = np.ravel(hdul['SCIDATA'].data[r,span:-span-1])
            dy = np.ravel(hdul['ERRDATA'].data[r,span:-span-1])
            q = np.ravel(hdul['QUALDATA'].data[r,span:-span-1])
        else:
            x = hdul['WAVEDATA_VAC_BARY'].data[row,span:-span-1]
            y = hdul['SCIDATA'].data[row,span:-span-1]
            dy = hdul['ERRDATA'].data[row,span:-span-1]
            q = hdul['QUALDATA'].data[row,span:-span-1]



        xmin, xmax = self._create_xmin_xmax(x)
        w = np.where(xmax-xmin > 0)
        x,xmin,xmax = x[w],xmin[w],xmax[w]
        y = y[w]/(xmax-xmin)#*10#au.nm/au.Angstrom
        dy = dy[w]/(xmax-xmin)#*10#au.nm/au.Angstrom
        q = q[w]
        argsort = np.argsort(x)
        x,xmin,xmax,y,dy,q = x[argsort],xmin[argsort],xmax[argsort],y[argsort],dy[argsort],q[argsort]

        w = np.where(q<1) #4**7)
        x,xmin,xmax,y,dy,q = x[w],xmin[w],xmax[w],y[w],dy[w],q[w]

        resol = []*len(x)
        xunit = au.Angstrom
        yunit = au.erg / au.Angstrom / au.cm**2 / au.s
        meta = hdr
        spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
        spec._t['quality'] = q
        return spec


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

    def _col_name(self, data, col):
        flag = '--'+col+'col'
        if self._gui._flags_cond(flag):
            col_name = self._gui._flags_extr(flag)
            if col_name in data.colnames:
                logging.info("I am using '%s' as '%s' column." % (col_name, col))
                return col_name
            else:
                logging.info("I can't find column '%s', I will use standard "
                             "column names instead." % col_name)
        if col == 'x': col_names = x_col_names
        if col == 'y': col_names = y_col_names
        if col == 'dy': col_names = dy_col_names
        if col == 'cont': col_names = cont_col_names
        col_sel = np.where([c in data.colnames for c in col_names])[0]
        try:
            col_sel = [col_sel[0]]
        except:
            pass
        if len(col_sel) == 1:
            return col_names[col_sel][0]
        else:
            return None


    def generic_spectrum(self, sess, hdul):
        """ Generic spectrum """
        logging.info(msg_format('generic'))
        hdr = hdul[0].header
        self._gui = sess._gui
        try:
            if len(hdul)>1:
                # MARZ spectra
                if len(hdul)==4 \
                    and hdul[1].header['EXTNAME']=='VARIANCE' \
                    and hdul[2].header['EXTNAME']=='WAVELENGTH' \
                    and hdul[3].header['EXTNAME']=='FIBRES':
                    x = hdul[2].data[0]
                    y = hdul[0].data[0]
                    dy = hdul[1].data[0]
                    cont = []
                    data = None
                else:
                    data_s = hdul[1].data
                    data = Table(data_s)
                    x_name = self._col_name(data, 'x')
                    y_name = self._col_name(data, 'y')
                    dy_name = self._col_name(data, 'dy')
                    cont_name = self._col_name(data, 'cont')
                    try:
                        x = np.ravel(data[x_name])
                        y = np.ravel(data[y_name])
                        dy = data[dy_name] if dy_name is not None \
                            else np.full(len(y), np.nan)
                        cont = data[cont_name] if cont_name is not None \
                            else np.full(len(y), np.nan)
                    except:
                        logging.error("I can't recognize columns.")
                        return 0

            else:
                data_s = hdul[0].data
                data = hdul[0].data
                x = data[0][:]
                y = data[1][:]
                dy = data[2][:]
                try:
                    cont = data[3][:]
                except:
                    cont = []

            # Import unit (if present)
            try:
                xunit = data_s.__dict__['_coldefs'][x_name]._unit
            except:
                xunit = None

            if xunit == None:
                xunit = au.nm
                if np.nanmax(x)>3000:
                    x = x*0.1
                if np.nanmax(x)<3:
                    x = x*1000
            try:
                yunit = au.Unit(data_s.__dict__['_coldefs'][y_name]._unit)
            except:
                yunit = au.erg/au.cm**2/au.s/au.Angstrom


            xmin, xmax = self._create_xmin_xmax(x)
            meta = hdr #{}

            # De-normalize
            norm_check = np.median(y)*np.max(y)
            if norm_check > 0.7 and norm_check < 1.3 and not all(np.isnan(cont)):
                y = y*cont
                dy = dy*cont

            spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta, cont=cont)
            if data is not None and hasattr(data, 'colnames'):
                for i,c in enumerate(data.colnames):
                    if c not in [x_name, y_name, dy_name, 'xmax', 'xmin']:
                        spec._t[c] = data[c]
                    #spec._t[c].unit = hdr1['TUNIT%i' % (i+1)]
            return spec
        except:
            return None


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


    def lrs_spectrum(self, hdul):
        """ TNG LRS spectrum """
        logging.info(msg_format('LRS'))
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
        meta = hdr
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
        y = data[:][0]#*data[:][3]
        #dy = data[:][1]#*data[:][3]
        dy = data[:][2]#*data[:][3]
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


    def xqr30_bosman(self, hdul):
        """ XQR-30 spectrum as formatted by Sarah Bosman """

        logging.info(msg_format('xqr30_bosman'))
        hdr = hdul[0].header

        data = Table(hdul[1].data)
        x = data['col1']
        y = data['col2']
        dy = data['col3']
        xunit = au.Angstrom
        yunit = au.dimensionless_unscaled
        xmin, xmax = self._create_xmin_xmax(x)
        meta = hdr #{}
        spec = Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

        #spec._t['cont'] = data['col5']
        spec._t['redside_PCA'] = data['col4']
        spec._t['blueside'] = data['col5']
        spec._t['blueside_1sigmal'] = data['col6']
        spec._t['blueside_1sigmau'] = data['col7']
        spec._t['blueside_2sigmal'] = data['col8']
        spec._t['blueside_2sigmau'] = data['col9']
        return spec


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

    def efosc2_spectrum(self, hdul):
        """ NTT EFOSC2 spectrum """
        logging.info(msg_format('NTT EFOSC2'))
        hdr = hdul[0].header
        data = hdul[0].data
        crval1 = hdr['CRVAL1']
        cdelt1 = hdr['CDELT1']
        naxis1 = hdr['NAXIS1']
        crpix1 = hdr['CRPIX1']
        y = data
        x = np.arange(crval1, crval1+naxis1*cdelt1, cdelt1)[:len(y)] - crpix1*cdelt1
        xmin, xmax = self._create_xmin_xmax(x)
        dy = np.full(len(y), np.nan)
        xunit = au.Angstrom
        yunit = au.erg/au.cm**2/au.s/au.Angstrom
        meta = hdr
        return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)

def sdss_spectrum(self, hdul):
    """SDSS spectrum"""
    logging.info(msg_format('SDSS'))

    def revIVar(x, m):
        if x == 0:
            return m
        return np.sqrt(1 / x)
    vectRevIVar = np.vectorize(revIVar)

    hdr = hdul[1].header
    data = np.array([np.array(i) for i in hdul[1].data])

    y = data[:, 0]
    x = (10 ** data[:, 1])
    xmin, xmax = self._create_xmin_xmax(x)
    dy = vectRevIVar(data[:, 2], max(y))

    xunit = au.Angstrom
    yunit = au.erg/au.cm**2/au.s/au.Angstrom
    meta = hdr
    return Spectrum(x, xmin, xmax, y, dy, xunit, yunit, meta)
