
import logging
import numpy as np
import unittest
from unittest.mock import MagicMock
from astrocook.core.session import SessionV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.recipes.features import RecipeFeaturesV2

class TestRecipeFeaturesV2(unittest.TestCase):
    def setUp(self):
        # Create a mock session
        self.session = MagicMock(spec=SessionV2)
        
        # Create a simple synthetic spectrum
        # Continuum = 1.0, Absorption feature at 500 nm with depth 0.5
        x = np.linspace(400, 600, 201) # 1 nm steps
        y = np.ones_like(x)
        # Add feature: Gaussian centered at 500 nm, sigma=2 nm, depth=0.5
        y -= 0.5 * np.exp(-0.5 * ((x - 500) / 2)**2)
        dy = np.ones_like(x) * 0.1
        
        # Mock spectrum object
        self.spec = MagicMock(spec=SpectrumV2)
        self.spec.x.value = x
        self.spec.y.value = y
        self.spec.dy.value = dy
        # basic attributes needed
        self.spec.has_aux_column.return_value = False
        
        self.session.spec = self.spec
        
        # Initialize recipe
        self.recipe = RecipeFeaturesV2(self.session)

    def test_compute_ew_linear_continuum(self):
        # Test basic EW calculation with linear continuum interpolation
        # Feature is roughly between 490 and 510
        # True EW of Gaussian = depth * sigma * sqrt(2*pi)
        # Depth=0.5, sigma=2 -> EW = 0.5 * 2 * 2.5066 = 2.5066 nm
        
        # Selection range: 490 to 510 (where y is ~1.0 at edges)
        xmin = 490.0
        xmax = 510.0
        
        result = self.recipe.compute_ew(xmin, xmax)
        
        self.assertIn('ew', result)
        ew = result['ew']
        centroid = result['centroid']
        
        print(f"Computed EW: {ew:.4f} nm")
        print(f"Computed Centroid: {centroid:.4f} nm")
        
        # Expected value approx 2.5 nm
        self.assertAlmostEqual(ew, 2.5066, delta=0.1) # 10% tolerance is generous for discrete integration
        self.assertAlmostEqual(centroid, 500.0, delta=0.5)

    def test_compute_ew_no_points(self):
        # Test handling of empty selection
        xmin = 700.0
        xmax = 710.0
        result = self.recipe.compute_ew(xmin, xmax)
        self.assertEqual(result, {})

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
