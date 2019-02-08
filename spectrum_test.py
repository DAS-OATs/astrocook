from astrocook.spectrum import Spectrum
from copy import deepcopy as dc
import unittest

x = [1, 2, 3]
xmin = [0.9, 1.8, 2.7]
xmax = [1.1, 2.2, 3.3]
y = [5, 6, 7]
dy = [0.3, 0.4, 0.5]
check_empty = []
spec_empty = Spectrum()
spec_full = Spectrum(x=x, xmin=xmin, xmax=xmax, y=y, dy=dy)

par_extract_region = {'xmin': 1.5, 'xmax': 2.5}
x = [2]
xmin = [1.8]
xmax = [2.2]
y = [6]
dy = [0.4]
check_extract_region = Spectrum(x=x, xmin=xmin, xmax=xmax, y=y, dy=dy)

class SpectrumTest(unittest.TestCase):

    def equal(self, spec, check):
        for s, c in zip(spec.t.colnames, check.t.colnames):
            self.assertListEqual(list(s), list(c))

    def test_extract_region(self):
        reg = spec_full.extract_region(**par_extract_region)
        self.equal(reg, check_extract_region)

if __name__ == '__main__':
    unittest.main()
