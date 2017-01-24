from astrocook import Line
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import testing as nput
from unittest import TestCase, TestLoader, TextTestRunner
import sys
    
class LineTest(TestCase):
    """Class for tests on line lists
    
    All tests on Line should be defined as method of this class.
    """
    
    def test___init___(self):
        """Test of Line.__init___
        
        A line list is initialized from a FITS file and plotted.
        """

        line = Line(file = 'FLINE_PRE.fits')
        x = line['bin']
        y = line['flux']
        plt.figure(1)
        plt.xlabel(x.unit)
        plt.ylabel(y.unit)
        plt.plot(x, y, 'o')
        plt.show()
        
    def test_misc(self):
        """Test of Line.misc
        
        A line list is initialized from a FITS file. The miscellaneous 
        information is updated with new entries and checked. 
        """
        
        line = Line(file = 'FLINE_PRE.fits')
        self.assertEqual(line.misc, {})        
        line.misc.update({'a': 1.})
        self.assertEqual(line.misc, {'a': 1.})        
        line.misc['b'] = 2
        self.assertEqual(line.misc, {'a': 1., 'b': 2})        
        line.misc.update(b = 3.)
        self.assertEqual(line.misc, {'a': 1., 'b': 3.})        
        
if __name__ == '__main__':
    """Run the tests"""

    suite = TestLoader().loadTestsFromTestCase(LineTest)
    TextTestRunner(verbosity = 2).run(suite)