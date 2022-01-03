from . import version
from .cookbook import Cookbook
from .defaults import Defaults
from .format import Format
from .functions import *
from .gui_log import GUILog
from .line_list import LineList
from .lmfit_model import load_model, save_model
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
from astropy.table import Column, Table
from collections import OrderedDict
from copy import deepcopy as dc
import json
from lmfit.model import Model
import logging
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.transforms as transforms
import numpy as np
import operator
import os
import pickle
#import dill as pickle
from scipy.signal import argrelmin
import shelve
import tarfile
import time


class Session(object):
    """ Class for sessions.

    A Session is a self-sufficient set of analysis operations."""

    def __init__(self,
                 gui=None,
                 path=None,
                 name=None,
                 spec=None,
                 spec_form=None,
                 nodes=None,
                 lines=None,
                 systs=None,
                 mods=None,
                 twin=False,
                 row=None,
                 slice=None):
        self._gui = gui
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
        self.log = GUILog(self._gui)
        self.defs = Defaults(self._gui)
        self._open_twin = twin
        self._row = row
        self._slice = slice
        self._clicks = []
        self._stats = False
        self._shade = False


    def _append(self, frame, append=True):
        if append and hasattr(self, frame.__name__):
            getattr(self, frame.__name__)._append(frame)
        else:
            setattr(self, frame.__name__, frame)

    def _data_iden(self, hdul, hdr):

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
            telesc = hdr['TELESCOP']
        except:
            telesc = 'undefined'

        try:
            hist = [i.split(' ') for i in str(hdr['HISTORY']).split('\n')]
            hist = [i for j in hist for i in j]
            if 'UVES_popler:' in hist and orig != 'Astrocook':
                instr = 'UVES'
                orig = 'POPLER'
            if 'XSHOOTER_REDUCE' in hist and orig != 'Astrocook':
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

        try:
            prefix = self.path.split('/')[-1][:2]
            #prefix = os.path.realpath(self.path).split('/')[-1][:2]
            if prefix == 'ql':
                orig = 'QUBRICS'
        except:
            pass

        if instr == None:
            logging.warning(msg_descr_miss('INSTRUME'))
        if catg == None:
            logging.warning(msg_descr_miss('HIERARCH ESO PRO CATG'))
        if orig == None:
            logging.warning(msg_descr_miss('ORIGIN'))

        logging.debug("Instrument: %s; origin: %s; category: %s."
                      % (instr, orig, catg))

        self._instr, self._catg, self._orig, self._telesc = instr, catg, orig, telesc

    def _other_open(self, hdul, hdr):
        format = Format()

        instr, catg, orig, telesc = self._instr, self._catg, self._orig, self._telesc

        # ESO ADP spectrum
        if orig == 'ESO' and hdr['ARCFILE'][:3]=='ADP':
            self.spec = format.eso_adp(hdul)

        # ESO-MIDAS spectrum
        if orig == 'ESO-MIDAS' and telesc != 'ESO-NTT':
            if len(hdul) == 1:
                #self.spec = format.eso_midas_image(hdul)
                self.spec = format.generic_spectrum(self, hdul)
            else:
                #self.spec = format.eso_midas_table(hdul)
                self.spec = format.generic_spectrum(self, hdul)

        # NTT spectra are parsed incorrectly, temp fix while Guido fixes the issue
        if orig == 'ESO-MIDAS' and telesc == 'ESO-NTT':
            self.spec = format.efosc2_spectrum(hdul)

        # ESPRESSO S1D spectrum
        if instr == 'ESPRESSO' and catg[0:3] == 'S1D':
            self.spec = format.espresso_s1d_spectrum(hdul)
            p = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/../'
            self.spec_form = format.espresso_spectrum_format(
                ascii.read(p+'espr_spec_form.dat'))

        # ESPRESSO S2D spectrum
        if instr == 'ESPRESSO' and catg[0:3] == 'S2D':
            self.spec = format.espresso_s2d_spectrum(hdul, row=self._row,
                                                     slice=self._slice)
            p = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/../'
            self.spec_form = format.espresso_spectrum_format(
                ascii.read(p+'espr_spec_form.dat'))
            if self._row is None:
                self._row = 0
            elif self._row < 169:
                self._row += 1
            else:
                self._row = None
            self._order = np.append(np.repeat(range(161,116,-1), 2),
                                    np.repeat(range(117,77,-1), 2))

        # ESPRESSO DAS spectrum
        if instr in ('ESPRESSO', 'UVES') and catg[1:5] == 'SPEC':
            self.spec = format.espresso_das_spectrum(hdul)
            p = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/../'
            self.spec_form = format.espresso_spectrum_format(
                ascii.read(p+'espr_spec_form.dat'))

        # FIRE spectrum
        if instr == 'FIRE':
            self.spec = format.firehose_spectrum(hdul)

        # FIRE spectrum
        if instr == 'MagE':
            self.spec = format.mage_spectrum(hdul)

        # QUBRICS spectrum
        if orig == 'QUBRICS':
            self.spec = format.qubrics_spectrum(hdul)

        # TNG LRS spectrum
        if instr == 'LRS' and orig == 'ESO-MIDAS':
            self.spec = format.lrs_spectrum(hdul)

        # UVES Spectrum
        if instr == 'UVES':
            if 'FLUXCAL_SCI' in self.path:
                hdul_err = fits.open(self.path.replace('FLUXCAL_SCI',
                                                       'FLUXCAL_ERRORBAR_SCI'))
                self.spec = format.uves_spectrum(hdul, hdul_err)

        # UVES POPLER spectrum
        if instr == 'UVES' and orig == 'POPLER':
            self.spec = format.uves_popler_spectrum(hdul)

        # WFCCD Spectrum
        if instr[:5] == 'WFCCD':
            self.spec = format.wfccd_spectrum(hdul)

        # XSHOOTER MERGE1D spectrum
        if instr == 'XSHOOTER' and 'MERGE1D' in catg.split('_'):
            if hdul[0].header['NAXIS'] == 0:
                self.spec = format.xshooter_vacbary_spectrum(hdul)
            else:
                self.spec = format.xshooter_merge1d_spectrum(hdul)

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

        # generic
        if instr == 'undefined' and orig == 'undefined' and catg == 'undefined':
            if self.path[-3:]=='txt' and len(Table(hdul[1].data).colnames)==9:
                self.spec = format.xqr30_bosman(hdul)
            else:
                self.spec = format.generic_spectrum(self, hdul)


    def open(self):

        dat = False
        if self.path[-3:] == 'acs':
            root = '/'.join(self.path.split('/')[:-1])
            #root =  '/'.join(os.path.realpath(self.path).split('/')[:-1])
            with tarfile.open(self.path) as arch:
                arch.extractall(path=root)
                try:
                    hdul = fits.open(self.path[:-4]+'_spec.fits')
                    hdr = hdul[1].header
                except:
                    dat = True
        elif self.path[-4:] == 'fits' or self.path[-7:] == 'fits.gz':
            hdul = fits.open(self.path)
            hdr = hdul[0].header
        else:
            t = Table(ascii.read(self.path))
            hdul = fits.HDUList([fits.PrimaryHDU(),
                                 fits.BinTableHDU.from_columns(np.array(t))])
            hdr = hdul[0].header

        self._data_iden(hdul, hdr)
        #print(self._orig, self._catg, self._instr)

        # Astrocook structures
        format = Format()
        only_constr = False
        fast = False
        if (self._orig[:9] == 'Astrocook' and self.path[-3:] == 'acs') or dat:
            for s in self.seq:
                try:
                    hdul = fits.open(self.path[:-4]+'_'+s+'.fits')
                    setattr(self, s, format.astrocook(hdul, s))
                    os.remove(self.path[:-4]+'_'+s+'.fits')
                    os.remove(self.path[:-4]+'_'+s+'.dat')
                except:
                    try:
                        data = ascii.read(self.path[:-4]+'_'+s+'.dat')
                        setattr(self, s, format.astrocook(data, s))
                        #os.remove(self.path[:-4]+'_'+s+'.dat')
                    except:
                        pass
                if s == 'systs':
                    try:
                        data = ascii.read(self.path[:-4]+'_'+s+'_mods.dat')
                    except:
                        data = None

                    if data is not None:
                        systs = getattr(self, 'systs')
                        data = ascii.read(self.path[:-4]+'_'+s+'_mods.dat')
                        os.remove(self.path[:-4]+'_'+s+'_mods.dat')
                        #print(data.info)
                        #data['id'] = object
                        #for i,id in enumerate(data['id']):
                        #    data['id'][i] = np.array(list(map(int, id[1:-1].split(','))))
                        #print(data.info)
                        setattr(systs, '_mods_t', data['z0', 'chi2r'])
                        systs._mods_t.remove_column('chi2r')

                        #systs._mods_t['mod'] = None
                        systs._mods_t['mod'] = np.empty(len(data), dtype=object)
                        systs._mods_t['chi2r'] = data['chi2r']
                        systs._mods_t['id'] = np.empty(len(data), dtype=object)
                        for i in range(len(data)):
                            #print(np.array(list(map(int, data['id'][i][1:-1].split(',')))))
                            #print(type(np.array(list(map(int, data['id'][i][1:-1].split(','))))))
                            #systs._mods_t['id'][i] = list(np.array(list(map(int, data['id'][i][1:-1].split(',')))))
                            systs._mods_t['id'][i] = list(map(int, data['id'][i][1:-1].split(',')))
                        funcdefs = {'convolve_simple': convolve_simple,
                                    'lines_voigt': lines_voigt,
                                    'psf_gauss': psf_gauss,
                                    'zero': zero}
                        for i,m in enum_tqdm(systs._mods_t, len(systs._mods_t),
                                             "session: Opening models"):
                        #for m in systs._mods_t:
                            name_mod_dat = self.path[:-4]+'_'+s+'_mods_%i.dat' % m['id'][0]
                            with open(name_mod_dat, 'rb') as f:
                                mod = pickle.load(f)
                            #setattr(mod.__init__, '_tmp', mod.func)

                            for attr in ['_lines', '_group', 'left', 'right']:
                                name_attr_dat = self.path[:-4]+'_'+s\
                                    +'_mods_%i_%s.dat' % (m['id'][0], attr)
                                setattr(mod, attr, load_model(name_attr_dat,
                                        funcdefs=funcdefs))
                                os.remove(name_attr_dat)
                            super(SystModel, mod).__init__(mod._group, Model(zero), operator.add)
                            #super(SystModel, mod).__init__(mod.left, mod.right, mod.op)
                            class_unmute(mod, Spectrum, self.spec)
                            m['mod'] = mod
                            os.remove(name_mod_dat)


                        for m in systs._mods_t['mod']:
                            for attr in ['_mods_t']:
                                setattr(m, attr, getattr(systs, attr))


                            #print(m.func)

                        only_constr = True
                        fast = True

            if self.spec is not None and self.systs is not None:
                self.cb._mods_recreate(only_constr=only_constr, fast=fast)
                self.cb._spec_update()
            try:
                os.remove(self.path[:-4]+'.json')
            except:
                pass

        else:
            self._other_open(hdul, hdr)

        #if self._gui._flags is not None \
        #    and '--systs' in [f[:7] for f in self._gui._flags]:
        if self._gui._flags_cond('--systs'):
            #paths = [f.split('=')[-1] for f in self._gui._flags if f[:7]=='--systs']
            #if len(paths)>1:
            #    logging.warning("You gave me too many system lists! I will "\
            #                    "load the first one.")
            path = self._gui._flags_extr('--systs')
            data = ascii.read(path)
            # Only Ly_a for now
            series = ['Ly_a']*len(data)
            func = ['voigt']*len(data)
            z = data['col1']
            dz = data['col3']
            logN = data['col4']
            dlogN = data['col5']
            b = data['col6']
            db = data['col7']
            self.cb.resol_est()
            resol = [np.nanmean(self.spec._t['resol'])]*len(data)
            chi2r = [np.nan]*len(data)
            id = range(len(data))
            out = SystList(func=func, series=series, z=z, dz=dz, logN=logN,
                           dlogN=dlogN, b=b, db=db, resol=resol, chi2r=chi2r,
                           id=id)
            out._t['z0'] = z
            self.spec._t['cont'] = 1
            setattr(self, 'systs', out)
            #print(out._t)
            self.cb._mods_recreate()
            self.cb._spec_update()


    def save(self, path):

        root = path[:-4]
        stem = pathlib.PurePath(path[:-4]).parts[-1]

        import warnings
        warnings.filterwarnings("ignore")

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
                    name_dat = root+'_'+s+'.dat'
                    obj = dc(getattr(self, s))
                    t = dc(obj._t)
                    if s == 'systs':
                        name_mods_dat = root+'_'+s+'_mods.dat'
                        mods_t = dc(obj._mods_t['z0', 'chi2r', 'id'])
                        ids = []
                        for id in obj._mods_t['id']:
                            id_list = [int(i) for i in id]
                            ids.append(json.dumps(id_list))
                        mods_t['id'] = ids


                    for c in t.colnames:
                        if type(t[c][0]) == np.int64:
                            pass
                        elif type(t[c][0]) == str or type(t[c][0]) == np.str_:
                            pass
                        elif type(t[c][0]) == OrderedDict:
                            pass
                        elif type(t[c][0]) == dict:
                            pass
                        else:
                            if c in ['logN', 'dlogN', 'b', 'db', 'resol', 'chi2r', \
                            'snr']:
                                format = '%3.3'
                            else:
                                format = '%3.7'
                            if np.abs(np.nanmedian(t[c]))<1e-7:# and t[c][0]!=0:
                                format += 'e'
                            else:
                                format += 'f'
                            t[c].info.format = format
                        #print(c, type(t[c][0]), np.abs(np.median(t[c])), format)

                    if s!='systs':
                        t['x'] = t['x'].to(au.nm)
                        t['xmin'] = t['xmin'].to(au.nm)
                        t['xmax'] = t['xmax'].to(au.nm)
                    del_list = []
                    for i, k in enumerate(obj._meta):
                        if k in forbidden_keywords or k[:5] in forbidden_keywords:
                            del_list.append(i)
                    for i in del_list[::-1]:
                        del obj._meta[i]
                    t.meta = dc(obj._meta)
                    #print(t.meta.comments)
                    t.meta['ORIGIN'] = 'Astrocook'
                    #t.meta['HIERARCH ASTROCOOK VERSION'] = version
                    #t.meta['HIERARCH ASTROCOOK STRUCT'] = s
                    if s == 'systs':
                        for i,(k,v) in enumerate(obj._constr.items()):
                            if v[0] in t['id']:
                                t.meta['HIERARCH AC CONSTR ID %i' % i] = v[0]
                                t.meta['HIERARCH AC CONSTR PAR %i' % i] = v[1]
                                t.meta['HIERARCH AC CONSTR VAL %i' % i] = v[2]
                    for c in t.colnames:
                        t[c].unit = au.dimensionless_unscaled
                    #print(t)
                    #t.write(name, format='fits', overwrite=True)
                    hdr = fits.Header(t.meta)
                    for c in t.meta:
                        try:
                            hdr.comments[c] = t.meta.comments[c]
                        except:
                            pass
                    phdu = fits.PrimaryHDU(header=hdr)
                    #print([Column(t[c]) for c in t.colnames])
                    #cols = fits.ColDefs([Column(c) for c in t.columns])
                    #cols = []
                    #for c in colnames:
                    #    cols.append(Columns)
                    thdu = fits.BinTableHDU(data=t, header=hdr)
                    hdul = fits.HDUList([phdu, thdu])
                    hdul.writeto(name, overwrite=True)
                    #print([t[c].format for c in t.colnames] )
                    try:
                        ascii.write(t, name_dat, names=t.colnames,
                        	        format='commented_header', overwrite=True)
                        arch.add(name, arcname=stem+'_'+s+'.fits')
                        arch.add(name_dat, arcname=stem+'_'+s+'.dat')
                        os.remove(name)
                        os.remove(name_dat)
                        logging.info("I've saved frame %s as %s."
                                     % (s, stem+'_'+s+'.fits/.dat'))
                    except:
                        logging.warning("I cannot save structure %s in ASCII "
                                        "format." % s)
                        arch.add(name, arcname=stem+'_'+s+'.fits')
                        os.remove(name)
                        logging.info("I've saved frame %s as %s."
                                     % (s, stem+'_'+s+'.fits'))
                    if s == 'systs':
                        ascii.write(mods_t, name_mods_dat,
                                    names=['z0', 'chi2r', 'id'],
                                    format='commented_header', overwrite=True)
                        arch.add(name_mods_dat, arcname=stem+'_'+s+'_mods.dat')
                        os.remove(name_mods_dat)

                        #name_mods_db = '%s.db' % (name_mods_dat[:-4])
                        #db = shelve.open(name_mods_db)

                        for i, r in enum_tqdm(obj._mods_t, len(obj._mods_t),
                                                "session: Saving models"):
                            id = r['id']
                            m = r['mod']
                            #print(m.__dict__)
                        #for id, m in obj._mods_t['id', 'mod']:
                            #find_class(m, Spectrum)
                            class_mute(m, Spectrum)
                            #print(m.func, m.left)

                            #mn = SystModel(m._spec, getattr(self, s))
                            #print(m.__dict__)
                            """
                            for attr in ['_series', '_vars', '_constr', '_z0',
                                         '_resol', '_defs', '_lines_pref',
                                         '_psf_pref', '_pars', '_xs',
                                         '_group_list', '_group_sel', '_ys',
                                         'op', '_prefix', '_param_root_names',
                                         'independent_vars', '_func_allargs',
                                         '_func_haskeywords', 'nan_policy',
                                         'opts', 'param_hints', '_param_names',
                                         'def_vals', '_name', '_xr', '_yr',
                                         '_wr', '_xf', '_yf', '_wf']:
                                setattr(mn, attr, getattr(m, attr))
                            """
                            for attr in ['_lines', '_group', 'left', 'right']:
                                name_attr_dat = '%s_%i_%s.dat' % (name_mods_dat[:-4], id[0], attr)
                                #if attr == 'func': print(getattr(m, attr).__dict__)
                                save_model(getattr(m, attr), name_attr_dat)
                                arch.add(name_attr_dat, arcname=stem+'_'+s+'_mods_%i_%s.dat' % (id[0], attr))
                                os.remove(name_attr_dat)

                            for attr in ['_mods_t']:
                                setattr(m, attr, attr)

                            for attr in ['_lines', '_group', 'left',
                                         'right', 'func']:
                                #setattr(mn, attr, None)
                                setattr(m, attr, None)
                            #print(mn.__dict__)
                            """
                            m._spec = ''
                            m.opts['spec'] = ''
                            m._group.opts['spec'] = ''
                            m._group.left.opts['spec'] = ''
                            m._group.left.right.opts['spec'] = ''
                            m._lines.opts['spec'] = ''
                            m._lines.right.opts['spec'] = ''
                            """
                            name_mod_dat = '%s_%i.dat' % (name_mods_dat[:-4], id[0])
                            #save_model(m, name_mod_dat)
                            with open(name_mod_dat, 'wb') as f:
                                #for a in m.__dict__:
                                pickle.dump(m, f, pickle.HIGHEST_PROTOCOL)

                            #db['mod'] = mn
                            arch.add(name_mod_dat, arcname=stem+'_'+s+'_mods_%i.dat' % id[0])
                            os.remove(name_mod_dat)
                        #arch.add(name_mods_db, arcname=stem+'_'+s+'_mods.db')
                        #os.remove(name_mods_db)
                        #db.close()


                    """
					arch.add(name, arcname=stem+'_'+s+'.fits')
                    arch.add(name_dat, arcname=stem+'_'+s+'.dat')

                    os.remove(name)
                    os.remove(name_dat)
                    logging.info("I've saved frame %s as %s."
                                 % (s, stem+'_'+s+'.fits'))
					"""
                #else:
                #    logging.warning("I haven't found any frame %s to save." % s)

            file = open(root+'.json', "w")
            n = file.write(self.log.str)
            file.close()
            arch.add(root+'.json', arcname=stem+'.json')
            os.remove(root+'.json')

    def save_pdf(self, path):

        """
             plt.figure(figsize=(3, 3))
             plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
             plt.title('Page One')
             plt.show()
             pdf.savefig()  # saves the current figure into a pdf page
             plt.close()
        """
        main = self._gui._graph_main
        graph = main._graph
        panel_n = 4
        length = 2
        x = self.spec.x.to(au.nm).value
        y = self.spec.y.value
        xmin = int(np.floor(np.nanmin(x/length))*length)
        xmax = int(np.ceil(np.nanmax(x/length))*length)
        xran = np.arange(xmin, xmax+length, length)
        page_n = int(np.ceil((xmax-xmin)/(length*panel_n)))
        i = 0
        with PdfPages(path) as pdf:
            #for p in range(page_n):
            for j, p in enum_tqdm(range(page_n), page_n-1,
                                  "session: Saving spectrum into PDF"):
                fig, axs = plt.subplots(panel_n,1, figsize=(11.69,8.27))
                fig.suptitle(self.name, size=18)
                for ax in axs:
                    try:
                        if (i+1)%panel_n == 0:
                            ax.text(0.48,-0.5, j+1, size=13, transform=ax.transAxes)

                        sel = np.where(np.logical_and(x>xran[i], x<xran[i+1]))
                        #ymin = np.floor(np.nanmin(y[sel]))
                        if main._norm:
                            yceil = 1
                            edge = 0.25
                        else:
                            yceil = np.ceil(np.nanmax(y[sel]))
                            edge = 0.1
                        ymin = -edge*yceil if yceil>0 else (1+edge)*yceil
                        ymax = (1+edge)*yceil if yceil>0 else -edge*yceil
                        #print(i,j,xran[i], xran[i+1],ymin, ymax)
                        ax.set_xlim(xran[i], xran[i+1])
                        ax.set_ylim(ymin, ymax)
                        graph._seq_core(self, main._norm, True, False, main, ax,
                                        cursor_list=['systs', 'cursor'])
                        if (i+1)%panel_n == 0:
                            ax.set_xlabel(self.spec.x.unit)
                        if main._norm:
                            ax.set_ylabel("Normalized flux")
                        else:
                            ax.set_ylabel(self.spec.y.unit)
                    except:
                        ax.set_frame_on(False)
                        ax.axes.get_xaxis().set_visible(False)
                        ax.axes.get_yaxis().set_visible(False)
                    i += 1

            #plt.show()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
