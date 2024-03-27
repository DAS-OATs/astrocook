from . import version
from .cookbook import Cookbook
from .defaults import Defaults
from .format import Format
from .functions import *
from .gui_log import GUILog
from .interv_list import IntervList
from .line_list import LineList
from .lmfit_model import load_model, save_model
from .message import *
#from .model import Model
from .spectrum import Spectrum
from .syst_list import SystList
from .syst_model import SystModel
from .feat_list import FeatList
#from .model_list import ModelList
from .vars import *
#from astropy import constants as ac
from astropy import units as au
from astropy.io import ascii, fits
from astropy.table import Column, Table
from collections import OrderedDict
from copy import deepcopy as dc
import glob
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
import shutil
import sys
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
                 intervs=None,
                 nodes=None,
                 lines=None,
                 systs=None,
                 mods=None,
                 twin=False,
                 row=None,
                 slice=None,
                 feats=None):
        self._gui = gui
        self.path = path
        self.name = name
        self.spec = spec
        self.spec_form = spec_form
        self.intervs = intervs
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
        self.feats = feats

        self._classes = {'spec': Spectrum, 'intervs': IntervList,
                         'lines': LineList, 'systs': SystList,
                         'mods': SystModel, 'feats': FeatList}


    def _append(self, frame, append=True):
        if append and hasattr(self, frame.__name__):
            getattr(self, frame.__name__)._append(frame)
        else:
            setattr(self, frame.__name__, frame)


    def _constr_from_mods(self):
        systs = self.systs
        systs._constr = {}
        for m in systs._mods_t['mod']:
            for p,v in m._pars.items():
                i, k = p.split('_')[-2:]
                i = int(i)
                if v.expr != None:
                    systs._constr[p] = (i, k, v.expr)
                if k in ['z', 'logN', 'b'] and not v.vary and v.expr==None:
                    systs._constr[p] = (i, k, None)


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

        # Impose generic format to Astrocook-processed data
        if orig=='Astrocook':
            catg = 'undefined'

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
            return 0

        # SDSS (should work for all release till 17)
        if telesc == "SDSS 2.5-M":
            self.spec = format.sdss_spectrum(hdul)
            return 0

        # ESO-MIDAS spectrum
        if orig == 'ESO-MIDAS' and telesc != 'ESO-NTT':
            if len(hdul) == 1:
                #self.spec = format.eso_midas_image(hdul)
                self.spec = format.generic_spectrum(self, hdul)
                return 0
            else:
                #self.spec = format.eso_midas_table(hdul)
                self.spec = format.generic_spectrum(self, hdul)
                return 0

        # NTT spectra are parsed incorrectly, temp fix while Guido fixes the issue
        if instr in ['EFOSC', 'LDSS3-'] or orig == 'ESO-MIDAS' and telesc == 'ESO-NTT':
            self.spec = format.efosc2_spectrum(hdul)
            return 0

        # ESPRESSO S1D spectrum
        if instr == 'ESPRESSO' and catg[0:3] == 'S1D':
            self.spec = format.espresso_s1d_spectrum(hdul)
            p = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/../'
            self.spec_form = format.espresso_spectrum_format(
                ascii.read(p+'espr_spec_form.dat'))
            return 0

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
            return 0

        # ESPRESSO DAS spectrum
        if instr in ('ESPRESSO', 'UVES') and catg[1:5] == 'SPEC':
            self.spec = format.espresso_das_spectrum(hdul)
            p = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + '/../'
            self.spec_form = format.espresso_spectrum_format(
                ascii.read(p+'espr_spec_form.dat'))
            return 0

        # FIRE spectrum
        if instr == 'FIRE':
            #self.spec = format.firehose_spectrum(hdul)
            self.spec = format.generic_spectrum(self, hdul)

        # FIRE spectrum
        if instr == 'MagE':
            self.spec = format.mage_spectrum(hdul)
            return 0

        # HARPN spectrum
        if instr == 'HARPN':
            self.spec = format.harpn_spectrum(hdul)
            return 0

        # QUBRICS spectrum
        if orig == 'QUBRICS':
            self.spec = format.qubrics_spectrum(hdul)
            return 0

        # TNG LRS spectrum
        if instr == 'LRS' and orig == 'ESO-MIDAS':
            self.spec = format.lrs_spectrum(hdul)
            return 0

        # UVES Spectrum
        if instr == 'UVES':
            if 'FLUXCAL_SCI' in self.path:
                hdul_err = fits.open(self.path.replace('FLUXCAL_SCI',
                                                       'FLUXCAL_ERRORBAR_SCI'))
                self.spec = format.uves_spectrum(hdul, hdul_err)
                return 0

        # UVES POPLER spectrum
        if instr == 'UVES' and orig == 'POPLER':
            self.spec = format.uves_popler_spectrum(hdul)
            return 0

        # WFCCD Spectrum
        if instr[:5] == 'WFCCD':
            self.spec = format.wfccd_spectrum(hdul)
            return 0

        # XSHOOTER MERGE1D spectrum
        if instr == 'XSHOOTER' and 'MERGE1D' in catg.split('_'):
            if hdul[0].header['NAXIS'] == 0:
                self.spec = format.xshooter_vacbary_spectrum(hdul)
                return 0
            else:
                self.spec = format.xshooter_merge1d_spectrum(hdul)
                return 0

        # XQR-30 spectrum
        if instr == 'XSHOOTER' and orig == 'XQR-30':
            self.spec = format.xqr30_spectrum(hdul, corr=self._open_twin)
            self._open_twin = not self._open_twin
            return 0

        # XSHOOTER DAS spectrum
        if instr == 'XSHOOTER' and catg[1:5] == 'SPEC':
            self.spec = format.xshooter_das_spectrum(hdul)
            return 0

        # XSHOOTER_REDUCE spectrum
        if instr == 'XSHOOTER' and orig == 'REDUCE':
            hdul_e = fits.open(self.path[:-5]+'e.fits')
            self.spec = format.xshooter_reduce_spectrum(hdul, hdul_e)
            return 0

        # generic
        if orig == 'undefined' and catg == 'undefined':
            if self.path[-3:]=='txt' and len(Table(hdul[1].data).colnames)==9:
                self.spec = format.xqr30_bosman(hdul)
                return 0
            else:
                self.spec = format.generic_spectrum(self, hdul)
                return 0

        # Astrocook spectrum-only
        if orig == 'Astrocook':
            self.spec = format.generic_spectrum(self, hdul)
            return 0


    def _rm_ac_temp(self, root):
        try:
            os.rmdir(root)
        except:
            logging.warning("I could not remove directory ac_temp/ (it was "\
                            "not empty).")


    def open(self):

        dat = False

        path = self.path
        parts = pathlib.PurePath(path[:-4]).parts
        stem = parts[-1]
        dir = parts[0].join(parts[0:-1])[1:]


        if self.path[-3:] == 'acs':
            root_super = '/'.join(self.path.split('/')[:-1])
            root = root_super+'/ac_temp/'
            with tarfile.open(self.path) as arch:
                #arch.extractall(path=root)
                arch.extractall(path=root)
                try:
                    try:
                        self._root_stem = root+stem
                        hdul = fits.open(self._root_stem+'_spec.fits')
                    except:
                        try:
                            g = glob.glob('ac_temp/*_spec.fits')[0]
                            s = root_super+'/'+g
                            hdul = fits.open(s)
                            self._root_stem = s[:-10]
                            logging.warning(
                                "The names of files within the .acs archive "\
                                "are different than expected: %s* instead of "\
                                "%s*." % (self._root_stem.split('/')[-1], stem))
                        except:
                            logging.error("I didn't find any *_spec.fits "
                                            "frame in the archive.")
                            self._rm_ac_temp(root)
                            return True

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

        # Astrocook structures
        format = Format()
        only_constr = False
        fast = False

        if (self._orig[:9] == 'Astrocook' and self.path[-3:] == 'acs') or dat:
            for s in self.seq:
                if s == 'feats':
                    try:
                        self._load(s, dir, stem, systs=self.systs)
                    except:
                        pass
                try:
                    hdul = fits.open(self._root_stem+'_'+s+'.fits')
                    setattr(self, s, format.astrocook(hdul, s))
                    os.remove(self._root_stem+'_'+s+'.fits')
                    os.remove(self._root_stem+'_'+s+'.dat')
                except:
                        try:
                            data = ascii.read(self._root_stem+'_'+s+'.dat')
                            setattr(self, s, format.astrocook(data, s))
                        except:
                            pass

                if s == 'systs':
                    try:
                        data = ascii.read(self._root_stem+'_'+s+'_mods.dat')
                    except:
                        data = None
                    if data is not None:
                        systs = getattr(self, 'systs')
                        data = ascii.read(self._root_stem+'_'+s+'_mods.dat')
                        os.remove(self._root_stem+'_'+s+'_mods.dat')
                        setattr(systs, '_mods_t', data['z0', 'chi2r'])
                        systs._mods_t.remove_column('chi2r')
                        systs._mods_t['mod'] = np.empty(len(data), dtype=object)
                        systs._mods_t['chi2r'] = data['chi2r']
                        systs._mods_t['id'] = np.empty(len(data), dtype=object)
                        for i in range(len(data)):
                            systs._mods_t['id'][i] = list(map(int, data['id'][i][1:-1].split(',')))

                        mods_t_ok = self._model_open(systs)
                        if mods_t_ok:
                            for m in systs._mods_t['mod']:
                                for attr in ['_mods_t']:
                                    setattr(m, attr, getattr(systs, attr))
                            self._constr_from_mods()
                            only_constr = True
                            fast = True
            if self.spec is not None and self.systs is not None:
                self.cb._mods_recreate(only_constr=only_constr, fast=fast)
                self.cb._spec_update()
                self.systs._dict_update(mods=True)

            if os.path.exists(self._root_stem+'.json'): os.remove(self._root_stem+'.json')

            self._rm_ac_temp(root)

        else:
            self._other_open(hdul, hdr)

        if self._gui._flags_cond('--systs'):
            path = self._gui._flags_extr('--systs')

            try:
                mode = self._gui._flags_extr('--mode')
            except:
                mode = 'std'

            # Only ascii for now
            logging.info("I'm using line list %s." % path)
            data = ascii.read(path)
            if mode == 'std':
                z = data['col1']
                dz = data['col3']
                logN = data['col4']
                dlogN = data['col5']
                b = data['col6']
                db = data['col7']
            if mode == 'viper':
                z = data['col1']/1215.67-1
                dz = [np.nan]*len(data)
                logN = data['col2']
                try:
                    dlogN = data['col5']
                except:
                    dlogN = [np.nan]*len(data)
                b = data['col3']
                try:
                    db = data['col6']
                except:
                    db = [np.nan]*len(data)

            # Only Ly_a for now
            series = ['Ly_a']*len(data)
            func = ['voigt']*len(data)
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
            self.cb._mods_recreate()
            self.cb._spec_update()



    def _model_open(self, systs):
        funcdefs = {'convolve_simple': convolve_simple,
                    'lines_voigt': lines_voigt,
                    'psf_gauss': self.spec.psf_gauss,
                    'zero': zero}

        mods_t_ok = True
        for i,m in enum_tqdm(systs._mods_t, len(systs._mods_t),
                             "session: Opening models"):
            try:
                name_mod_dat = self._root_stem+'_systs_mods_%i.dat' % m['id'][0]
                with open(name_mod_dat, 'rb') as f:
                    mod = pickle.load(f)

                for attr in ['_lines', '_group', 'left', 'right']:
                    name_attr_dat = self._root_stem+'_systs_mods_%i_%s.dat' % (m['id'][0], attr)
                    setattr(mod, attr, load_model(name_attr_dat,
                            funcdefs=funcdefs))
                    os.remove(name_attr_dat)
                super(SystModel, mod).__init__(mod._group, Model(zero), operator.add)
                class_unmute(mod, Spectrum, self.spec)
                m['mod'] = mod
                os.remove(name_mod_dat)
            except:
                mods_t_ok = False
        return mods_t_ok


    def _load(self, struct, dir, stem, **kwargs):

        new_dir = dir+'/'+stem+'_'+struct+'/'

        s = self._classes[struct]()
        """
        for file in os.listdir(new_dir):
            with open(new_dir+file, 'rb') as f:
                feats._l.append(pickle.load(f))
        """
        s._load(new_dir, **kwargs)
        setattr(self, struct, s)
        shutil.rmtree(new_dir, ignore_errors=True)
        logging.info("I loaded %s from %s.acs." % (struct, stem))


    def _save(self, struct, dir, stem, arch):
        if not hasattr(self, struct) or getattr(self, struct) is None:
            return None

        new_dir = dir+'/'+stem+'_'+struct+'/'
        try:
            shutil.rmtree(new_dir, ignore_errors=True)
            os.mkdir(new_dir)
        except:
            os.mkdir(new_dir)

        """
        l = self.feats._l

        for i, o in enumerate(l):
            with open(new_dir+'%04i.dat' % i, 'wb') as f:
                #for a in m.__dict__:
                pickle.dump(o, f, pickle.HIGHEST_PROTOCOL)
        """
        getattr(self, struct)._save(new_dir)

        arch.add(new_dir, arcname=stem+'_'+struct+'/')
        shutil.rmtree(new_dir, ignore_errors=True)
        logging.info("I've saved %s in %s.acs." % (struct, stem))


    def _save(self, struct, dir, stem, arch):
        if not hasattr(self, struct) or getattr(self, struct) is None:
            return None

        if dir!='':
            new_dir = dir+'/'+stem+'_'+struct+'/'
        else:
            new_dir = stem+'_'+struct+'/'
        try:
            shutil.rmtree(new_dir, ignore_errors=True)
            os.mkdir(new_dir)
        except:
            os.mkdir(new_dir)

        """
        l = self.feats._l

        for i, o in enumerate(l):
            with open(new_dir+'%04i.dat' % i, 'wb') as f:
                #for a in m.__dict__:
                pickle.dump(o, f, pickle.HIGHEST_PROTOCOL)
        """
        getattr(self, struct)._save(new_dir)

        arch.add(new_dir, arcname=stem+'_'+struct+'/')
        shutil.rmtree(new_dir, ignore_errors=True)
        logging.info("I've saved %s in %s.acs." % (struct, stem))


    def save(self, path, models=False):

        root = path[:-4]
        parts = pathlib.PurePath(path[:-4]).parts
        stem = parts[-1]
        dir = parts[0].join(parts[0:-1])[1:]

        import warnings
        warnings.filterwarnings("ignore")

        with tarfile.open(root+'.acs', 'w:gz') as arch:
            for s in self.seq:
                if s is 'feats':
                    self._save(s, dir, stem, arch)
                elif hasattr(self, s) and getattr(self, s) is not None:
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
                    #print(getattr(self, s).__dict__)
                    """
                    try:
                        mods = getattr(self, s)._mods_t
                        w = []
                        for m in mods:
                            w.append(744 in m['id'])

                        print(w)
                        mod = mods[w]
                        print(mod['id'][0])
                        print(mod['mod'][0].__dict__)
                        mod['mod'][0]._pars.pretty_print()
                    except:
                        pass
                    """
                    try:
                        obj = dc(getattr(self, s))
                    except:
                        obj = getattr(self, s)
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
                        try:
                            t['x'] = t['x'].to(au.nm)
                            t['xmin'] = t['xmin'].to(au.nm)
                            t['xmax'] = t['xmax'].to(au.nm)
                        except:
                            t['x'] = t['x'].to(au.km/au.s)
                            t['xmin'] = t['xmin'].to(au.km/au.s)
                            t['xmax'] = t['xmax'].to(au.km/au.s)
                    del_list = []
                    if hasattr(obj, '_meta'):
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
                    else:
                        hdr = fits.Header()

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
                        if not models: break

                        fail = []
                        for i, r in enum_tqdm(obj._mods_t, len(obj._mods_t),
                                                "session: Saving models"):
                            id = r['id']
                            m = r['mod']
                            try:
                                sys.setrecursionlimit(10000)
                                try:
                                    class_mute(m, Spectrum)
                                except:
                                    self.cb._mods_recreate(verbose=False)
                                    class_mute(m, Spectrum)
                                for attr in ['_lines', '_group', 'left', 'right']:
                                    name_attr_dat = '%s_%i_%s.dat' % (name_mods_dat[:-4], id[0], attr)
                                    save_model(getattr(m, attr), name_attr_dat)
                                    arch.add(name_attr_dat, arcname=stem+'_'+s+'_mods_%i_%s.dat' % (id[0], attr))
                                    os.remove(name_attr_dat)

                                for attr in ['_mods_t']:
                                    setattr(m, attr, attr)

                                attr_save = {}
                                for attr in ['_lines', '_group', 'left',
                                             'right', 'func']:
                                    attr_save[attr] = dc(getattr(m, attr))
                                    setattr(m, attr, None)
                                name_mod_dat = '%s_%i.dat' % (name_mods_dat[:-4], id[0])
                                with open(name_mod_dat, 'wb') as f:
                                    #for a in m.__dict__:
                                    pickle.dump(m, f, pickle.HIGHEST_PROTOCOL)

                                for attr in ['_lines', '_group', 'left',
                                             'right', 'func']:
                                    setattr(m, attr, attr_save[attr])
                                class_unmute(m, Spectrum, self.spec)

                                arch.add(name_mod_dat, arcname=stem+'_'+s+'_mods_%i.dat' % id[0])
                                os.remove(name_mod_dat)
                                sys.setrecursionlimit(1000)

                            except:
                                fail.append(id)

                        if fail != []:
                            logging.warning("I could not serialize %i out of %i "
                                            "models. They were not saved:" \
                                            % (len(fail), len(obj._mods_t)))
                            logging.warning("Here's the list of system IDs for "
                                            "the models that weren't saved:")
                            logging.warning(fail)

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
