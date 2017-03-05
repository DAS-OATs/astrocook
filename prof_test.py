from astrocook import ProfGauss
from astropy import units as u
import matplotlib.pyplot as plt
from numpy import testing as nput
from unittest import TestCase, TestLoader, TextTestRunner
import sys
    
class ProfGaussTest(TestCase):
    """Class for tests on profiles
    
    All tests on ProfGauss should be defined as method of this class.
    """
    
    def test___init___(self):
        """Test of ProfGauss.__init___
        
        A gaussian profile is created with default parameters and plotted.
        """

        prof = ProfGauss()
        x = prof['bin']
        y = prof['flux']
        plt.figure(1)
        plt.xlabel(x.unit)
        plt.ylabel(y.unit)
        plt.plot(x, y)
        plt.show()
        
if __name__ == '__main__':
    """Run the tests"""

    suite = TestLoader().loadTestsFromTestCase(ProfGaussTest)
    TextTestRunner(verbosity = 2).run(suite)