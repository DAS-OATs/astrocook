from . import Cont, Line, Spec1DReader, System
from astropy.io import fits
from astropy.table import Table
import difflib
import numpy as np
import os
import tarfile

class IO():
    def __init__(self):
        """ Constructor for the Spec1DReader class. """
        pass

    def acs_read(self, name, path):
        """ Read an astrocook session from a tar.gz archive """

        acs = IO()
        
        with tarfile.open(name) as arch:
            arch.extractall(path=path)
            
            # Spectrum (required)
            acs.spec_name = name[:-4]+'_spec.fits'
            #acs.spec_name = 'spec.fits'
            print acs.spec_name
            acs.spec = acs.spec_read(acs.spec_name)
            os.remove(acs.spec_name)
            
            # Line list
            try:
                acs.line_name = name[:-4]+'_line.fits'
                #acs.line_name = 'line.fits'
                acs.line = acs.line_read(acs.line_name)
                os.remove(acs.line_name)
                os.remove(acs.line_name[:-10]+'_mins.fits')
                os.remove(acs.line_name[:-10]+'_maxs.fits')
                os.remove(acs.line_name[:-10]+'_exts.fits')
            except:
                pass

            # Continuum
            try:
                acs.cont_name = name[:-4]+'_cont.fits'
                acs.cont = acs.cont_read(acs.cont_name)
                os.remove(acs.cont_name)
            except:
                pass


            # System list
            try:
                acs.syst_name = name[:-4]+'_syst.fits'
                acs.syst = acs.syst_read(acs.syst_name)
                os.remove(acs.syst_name)
            except:
                pass

        return acs

    def acs_write(self, acs, name, path, overwrite=True):
        """ Write an astrocook session onto a tar.gz archive """
        ndiff = difflib.ndiff(name, path)
        diff = ""
        for d in ndiff:
            if d[0] == '-':
                diff += d[-1]
        with tarfile.open(name, 'w:gz') as arch:
            spec_name = name[:-4] + '_spec.fits'
            spec_arcname = diff[:-4] + '_spec.fits'
            #self.spec_write(acs.spec, spec_name)
            acs.spec.save(spec_name)
            arch.add(spec_name, arcname=spec_arcname)
            os.remove(spec_name)

            try:
                line_name = name[:-4] + '_line.fits'
                line_arcname = diff[:-4] + '_line.fits'
                self.line_write(acs.line, line_name)
                arch.add(line_name, arcname=line_arcname)
                arch.add(line_name[:-10]+'_mins.fits',
                         arcname=line_arcname[:-10]+'_mins.fits')
                arch.add(line_name[:-10]+'_maxs.fits',
                         arcname=line_arcname[:-10]+'_maxs.fits')
                arch.add(line_name[:-10]+'_exts.fits',
                         arcname=line_arcname[:-10]+'_exts.fits')
                os.remove(line_name)
                os.remove(line_name[:-10]+'_mins.fits')
                os.remove(line_name[:-10]+'_maxs.fits')
                os.remove(line_name[:-10]+'_exts.fits')
            except:
                pass

            try:
                cont_name = name[:-4] + '_cont.fits'
                cont_arcname = diff[:-4] + '_cont.fits'
                self.cont_write(acs.cont, cont_name)
                arch.add(cont_name, arcname=cont_arcname)
                os.remove(cont_name)
            except:
                pass

            try:
                syst_name = name[:-4] + '_syst.fits'
                syst_arcname = diff[:-4] + '_syst.fits'
                self.syst_write(acs.syst, syst_name)
                arch.add(syst_name, arcname=syst_arcname)
                os.remove(syst_name)
            except:
                pass

                
    def cont_read(self, name):
        """ Read an astrocook continuum model from a FITS frame """
        
        hdul = fits.open(name)
        data = hdul[1].data
        names = np.array(hdul[1].columns.names)
        units = np.array(hdul[1].columns.units)
        xunit = units[np.where(names == 'X')][0]
        yunit = units[np.where(names == 'Y')][0]
        x = data['X']
        y = data['Y']            
        dy = data['DY']
        
        cont = Cont(x=x, y=y, dy=dy, xunit=xunit, yunit=yunit)

        return cont
        
    def cont_write(self, cont, name, overwrite=True):
        """ Write an astrocook continuum model onto a FITS frame """

        cont.t.write(name, format='fits', overwrite=overwrite)
    
    def line_read(self, name):
        """ Read an astrocook line list from a FITS frame """
        
        hdul = fits.open(name)
        data = hdul[1].data
        names = np.array(hdul[1].columns.names)
        units = np.array(hdul[1].columns.units)
        xunit = units[np.where(names == 'X')][0]
        yunit = units[np.where(names == 'Y')][0]
        x = data['X']
        xmin = data['XMIN']
        xmax = data['XMAX']
        y = data['Y']            
        dy = data['DY']
        
        line = Line(x=x, xmin=xmin, xmax=xmax, y=y, dy=dy, xunit=xunit,
                    yunit=yunit)

        # Extrema
        try:
            line._mins = Table.read(name[:-10]+'_mins.fits')
            line._maxs = Table.read(name[:-10]+'_maxs.fits')
            line._exts = Table.read(name[:-10]+'_exts.fits')
        except:
            pass
        
        return line

    def line_write(self, line, name, overwrite=True):
        """ Write an astrocook line list onto a FITS frame """

        line.t.write(name, format='fits', overwrite=overwrite)

        # Table with extrema
        line._mins.write(name[:-10]+'_mins.fits')
        line._maxs.write(name[:-10]+'_maxs.fits')
        line._exts.write(name[:-10]+'_exts.fits')

    def spec_read(self, name):
        """ Read a spectrum from a FITS frame """
        try:
            spec = Spec1DReader().ac(name)
        except:
            spec = Spec1DReader().xq(name)

        return spec
            
    def spec_write(self, spec, name, overwrite=True):
        """ Write an astrocook spectrum onto a FITS frame """

        spec.t.write(name, format='fits', overwrite=overwrite)
    
    def syst_read(self, name):
        """ Read an astrocook system from a FITS frame """

        hdul =  fits.open(syst_name)
        data = hdul[1].data
        names = np.array(hdul[1].columns.names)
        units = np.array(hdul[1].columns.units)
        Nunit = units[np.where(names == 'N')][0]
        bunit = units[np.where(names == 'B')][0]
        series = data['SERIES']
        z = data['Z']
        N = data['N']
        b = data['B']
        btur = data['BTUR']
        dz = data['DZ']
        dN = data['DN']
        db = data['DB']
        dbtur = data['DBTUR']
        vary = data['VARY']
        expr = data['EXPR']
        self.syst = System(
            spec=self.spec, cont=self.cont,
            series=series,
            z=z, N=N, b=b, btur=btur,
            dz=dz, dN=dN, db=db, dbtur=dbtur,
            vary=vary, expr=expr, Nunit=Nunit, bunit=bunit)
        syst.t.write(name, format='fits', overwrite=overwrite)

    def syst_write(self, syst, name, overwrite=True):
        """ Write an astrocook line list onto a FITS frame """

        syst.t.write(name, format='fits', overwrite=overwrite)
        
