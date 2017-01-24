from astrocook import Spec
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import testing as nput
from unittest import TestCase, TestLoader, TextTestRunner
import sys
    
class SpecTest(TestCase):
    """Class for tests on spectra
    
    All tests on spec.py should be defined as method of this class.
    """
    
    def test___init___(self):
        """Test of Spec.__init___
        
        A spectrum is initialized from a FITS file and plotted.
        """

        spec = Spec(file = 'RSPEC_1.fits')
        x = spec['bin']
        y = spec['flux']
        plt.figure(1)
        plt.xlabel(x.unit)
        plt.ylabel(y.unit)
        plt.plot(x, y)
        plt.show()
        
    def test_nexp(self):
        """Test of Spec.nexp
        
        A spectrum is initialized from a FITS file. The number of exposures is 
        checked.
        """
        
        spec = Spec(file = 'RSPEC_1.fits')
        self.assertEqual(spec.nexp, 1)        
        
    def test_exptime(self):
        """Test of Spec.exptime
        
        A spectrum is initialized from a FITS file. The exposure times are 
        checked, updated and checked again.
        """
        
        spec = Spec(file = 'RSPEC_1.fits')
        self.assertEqual(spec.exptime, [])        
        spec.exptime = [900.0 * u.s]
        self.assertEqual(spec.exptime, [900.0 * u.s])        
        spec.exptime = [900.0 * u.s, 600.0 * u.s]
        self.assertEqual(spec.exptime, [900.0 * u.s])        
        
    def test_misc(self):
        """Test of Spec.misc
        
        A spectrum is initialized from a FITS file. The miscellaneous 
        information is updated with new entries and checked. 
        """
        
        spec = Spec(file = 'RSPEC_1.fits')
        self.assertEqual(spec.misc, {})        
        spec.misc.update({'a': 1.})
        self.assertEqual(spec.misc, {'a': 1.})        
        spec.misc['b'] = 2
        self.assertEqual(spec.misc, {'a': 1., 'b': 2})        
        spec.misc.update(b = 3.)
        self.assertEqual(spec.misc, {'a': 1., 'b': 3.})        
        
if __name__ == '__main__':
    """Run the tests"""

    suite = TestLoader().loadTestsFromTestCase(SpecTest)
    TextTestRunner(verbosity = 2).run(suite)