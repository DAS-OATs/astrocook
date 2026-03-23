# import pytest
import numpy as np
import astropy.units as au
import logging

# Import V2 Core
from astrocook.core.structures import (
    SpectrumDataV2, DataColumnV2, ComponentDataV2, SystemListDataV2
)
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2

def create_synthetic_spectrum():
    """Generates a spectrum with a CIV doublet and some noise."""
    z_true = 2.0
    lambda_0 = 1548.2040
    lambda_obs_center = lambda_0 * (1 + z_true)
    
    # Create 500 pixels covering +/- 10 Angstroms
    x_arr = np.linspace(lambda_obs_center - 10, lambda_obs_center + 10, 500)
    
    # 2. Create True Model (Physics)
    true_comp = ComponentDataV2(
        id=0, z=z_true, dz=0, logN=14.0, dlogN=0, b=25.0, db=0, 
        series='CIV' 
    )
    
    dummy_list_data = SystemListDataV2(components=[true_comp])
    dummy_sys = SystemListV2(dummy_list_data)
    
    x_col = DataColumnV2(x_arr, au.Angstrom)
    y_dummy = DataColumnV2(np.ones_like(x_arr), au.dimensionless_unscaled)
    dy_dummy = DataColumnV2(np.ones_like(x_arr)*0.1, au.dimensionless_unscaled)
    spec_data = SpectrumDataV2(x=x_col, xmin=x_col, xmax=x_col, y=y_dummy, dy=dy_dummy)
    spec_v2 = SpectrumV2(spec_data)
    
    fitter = VoigtFitterV2(spec_v2, dummy_sys)
    p_true = np.array([z_true, 14.0, 25.0]) # z, logN, b (btur is frozen)
    true_flux = fitter._compute_model(p_true)
    
    # 3. Add Noise
    noise_level = 0.05
    noise = np.random.normal(0, noise_level, size=len(x_arr))
    y_obs = true_flux + noise
    dy_obs = np.ones_like(y_obs) * noise_level
    
    # 4. Return SpectrumV2 Object
    y_col = DataColumnV2(y_obs, au.dimensionless_unscaled)
    dy_col = DataColumnV2(dy_obs, au.dimensionless_unscaled)
    aux = {'cont': DataColumnV2(np.ones_like(x_arr), au.dimensionless_unscaled)}
    
    final_data = SpectrumDataV2(
        x=x_col, xmin=x_col, xmax=x_col, y=y_col, dy=dy_col, aux_cols=aux
    )
    return SpectrumV2(final_data)

def test_fitting_workflow():
    """
    Smoke test for the Frequentist Fitter (VoigtFitterV2).
    Verifies parameter recovery and compute_model_flux after refactor.
    """
    # 1. Get Data
    spectrum = create_synthetic_spectrum()
    
    # 2. Initial Guess
    z_true = 2.0
    logN_true = 14.0
    b_true = 25.0
    
    guess_comp = ComponentDataV2(
        id=1, 
        z=z_true + 0.0005,
        dz=None,
        logN=13.5, # Default guess
        dlogN=None,
        b=10.0,    # Default guess (required to trigger smart_guess)
        db=None,
        series='CIV'
    )
    
    syst_list = SystemListV2(SystemListDataV2(components=[guess_comp]))
    fitter = VoigtFitterV2(spectrum, syst_list)
    
    # 3. Run Fit
    new_syst_list, final_model_flux, result = fitter.fit(z_window_kms=100.0, verbose=0)
    
    # 4. Assertions
    assert result.success, "Fitter failed to converge"
    
    fitted = new_syst_list.components[0]
    
    # Tolerances (considering SNR=20)
    assert abs(fitted.z - z_true) < 1e-4, f"z recovery failed: {fitted.z} vs {z_true}"
    assert abs(fitted.logN - logN_true) < 0.2, f"logN recovery failed: {fitted.logN} vs {logN_true}"
    assert abs(fitted.b - b_true) < 5.0, f"b recovery failed: {fitted.b} vs {b_true}"

    # Verify model flux calculation
    x_calc, y_calc = fitter.compute_model_flux()
    assert len(y_calc) == len(spectrum.x.value)
    assert np.all(y_calc >= 0)

def test_adding_component_to_existing_system():
    """
    Reproduces the scenario where a user adds a new component to an existing fit.
    This tests if the fitter correctly handles index mapping when p_free is a 
    subset of the full parameters (i.e. partial fit).
    """
    spectrum = create_synthetic_spectrum()
    
    # 1. Existing System (Frozen/Inactive)
    comp1 = ComponentDataV2(
        id=1, z=1.5, dz=0, logN=13.0, dlogN=0, b=20.0, db=0, series='CIV'
    )
    
    # 2. New Component (Active/Free)
    comp2 = ComponentDataV2(
        id=2, z=2.0, dz=None, logN=13.5, dlogN=None, b=10.0, db=None, series='CIV'
    )
    
    syst_list = SystemListV2(SystemListDataV2(components=[comp1, comp2]))
    
    # Simulate GUI selection: only comp2 is active
    syst_list.constraint_model._active_uuids = [comp2.uuid]
    # Update constraints to freeze comp1
    # We must explicitly freeze comp1 as the GUI would do.
    # The default SystemListV2 constructor leaves dz=None components as free.
    # Here, comp1 has dz=0, which our constraints interpret as frozen.
    
    fitter = VoigtFitterV2(spectrum, syst_list)
    
    try:
        new_syst_list, final_model, res = fitter.fit(z_window_kms=100.0, verbose=0)
        assert res.success
    except IndexError as e:
        import pytest
        pytest.fail(f"IndexError occurred during partial fit: {e}")

if __name__ == "__main__":
    test_fitting_workflow()
    test_adding_component_to_existing_system()