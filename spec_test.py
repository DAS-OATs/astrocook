from astrocook import Spec
from astropy import units as u
from astropy.io import fits
from numpy import testing as nput
from unittest import TestCase, TestLoader, TextTestRunner
import sys
    
class SpecTest(TestCase):
    """Class for tests on spectra"""
    
    file = 'RSPEC.fits'

    def test___init___empty(self):
        spec = Spec()

    def test___init___from_fits(self):
        spec = Spec(file = self.file)

    def test_to_redshift(self):
    
        z = 3.1

        # Spectrum
        spec = Spec(file = self.file)
        data = fits.open(self.file)[1].data
        check_bin = data.field('WAVEL') * (1.0 + z)
        check_size = data.field('PIXSIZE') * (1.0 + z)
        spec.to_redshift(z)
        nput.assert_array_equal(spec['bin'], check_bin * u.nm)
        nput.assert_array_equal(spec['size'], check_size * u.nm)
        
if __name__ == '__main__':
    """Run the tests"""

    suite = TestLoader().loadTestsFromTestCase(SpecTest)
    TextTestRunner(verbosity = 2).run(suite)