import numpy as np
import astropy.units as au
import logging
import matplotlib.pyplot as plt

# Import V2 Core
from astrocook.core.structures import (
    SpectrumData, DataColumn, ComponentData, SystemListData
)
from astrocook.core.spectrum import Spectrum
from astrocook.core.system_list import SystemList
from astrocook.fitting.voigt_fitter import VoigtFitter

# Configure Logging
logging.basicConfig(level=logging.INFO)

def create_synthetic_spectrum():
    """Generates a spectrum with a CIV doublet and some noise."""
    print("--- Generating Synthetic Data ---")
    
    # 1. Create Grid centered on CIV 1548 at z=2.0
    z_true = 2.0
    lambda_0 = 1548.2040
    lambda_obs_center = lambda_0 * (1 + z_true)
    
    # Create 500 pixels covering +/- 10 Angstroms
    x_arr = np.linspace(lambda_obs_center - 10, lambda_obs_center + 10, 500)
    
    # 2. Create True Model (Physics)
    # We cheat and use the Fitter's own physics engine to generate the truth
    # True Params: z=2.0, logN=14.0, b=25.0
    true_comp = ComponentData(
        id=0, z=z_true, dz=0, logN=14.0, dlogN=0, b=25.0, db=0, 
        series='CIV' # This triggers both 1548 and 1550
    )
    
    # Create a dummy system list for generation
    dummy_list_data = SystemListData(components=[true_comp])
    dummy_sys = SystemList(dummy_list_data)
    
    # Create a dummy spectrum container
    x_col = DataColumn(x_arr, au.Angstrom)
    y_dummy = DataColumn(np.ones_like(x_arr), au.dimensionless_unscaled)
    dy_dummy = DataColumn(np.ones_like(x_arr)*0.1, au.dimensionless_unscaled)
    spec_data = SpectrumData(x=x_col, xmin=x_col, xmax=x_col, y=y_dummy, dy=dy_dummy)
    spec_v2 = Spectrum(spec_data)
    
    # Use Fitter to generate true flux
    fitter = VoigtFitter(spec_v2, dummy_sys)
    # We need to hack the constraint model to give us the P vector for the true component
    # Since we are just generating, we can calculate manually or use the internal method if accessible
    # Let's trust the fitter logic:
    # We create a fitter, but we need to generate the model from the *parameters*, not fit.
    # VoigtFitter._compute_model uses the constraints vector.
    
    _, true_flux = fitter.compute_model_flux()
    
    # 3. Add Noise (fixed seed: the test must be deterministic)
    noise_level = 0.05
    rng = np.random.default_rng(42)
    noise = rng.normal(0, noise_level, size=len(x_arr))
    y_obs = true_flux + noise
    dy_obs = np.ones_like(y_obs) * noise_level
    
    # 4. Return Spectrum Object
    y_col = DataColumn(y_obs, au.dimensionless_unscaled)
    dy_col = DataColumn(dy_obs, au.dimensionless_unscaled)
    
    # Add a continuum column (all ones)
    aux = {'cont': DataColumn(np.ones_like(x_arr), au.dimensionless_unscaled)}
    
    final_data = SpectrumData(
        x=x_col, xmin=x_col, xmax=x_col, y=y_col, dy=dy_col, aux_cols=aux
    )
    return Spectrum(final_data)

def _run_fit():
    
    # 1. Get Data
    spectrum = create_synthetic_spectrum()
    
    # 2. Create a "Bad" Initial Guess (User clicked slightly wrong)
    # True was z=2.0, logN=14.0, b=25.0
    # Guess is z=2.0001, logN=13.5, b=40.0
    print("\n--- Initializing System List with Guess ---")
    guess_comp = ComponentData(
        id=1, 
        z=2.0001,  # Slightly shifted (within the default z window of 20 km/s)
        dz=None,
        logN=13.5, # Too weak
        dlogN=None,
        b=40.0,    # Too broad
        db=None,
        series='CIV'
    )
    
    syst_list = SystemList(SystemListData(components=[guess_comp]))
    
    print(f"Guess: z={guess_comp.z:.5f}, logN={guess_comp.logN:.2f}, b={guess_comp.b:.2f}")
    
    # 3. Initialize Fitter
    print("\n--- Initializing Fitter ---")
    fitter = VoigtFitter(spectrum, syst_list)
    
    # 4. Run Fit
    print("--- Running Fit (Levenberg-Marquardt) ---")
    new_syst_list, _, result = fitter.fit(verbose=2)
    
    # 5. Inspect Results
    fitted_comp = new_syst_list.components[0]
    
    print("\n--- Fit Results ---")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Cost: {result.cost:.4f}")
    
    print("-" * 30)
    print(f"{'Param':<10} | {'True':<10} | {'Guess':<10} | {'Fitted':<10}")
    print("-" * 30)
    print(f"{'z':<10} | {2.0:<10.5f} | {guess_comp.z:<10.5f} | {fitted_comp.z:<10.5f}")
    print(f"{'logN':<10} | {14.0:<10.2f} | {guess_comp.logN:<10.2f} | {fitted_comp.logN:<10.2f}")
    print(f"{'b':<10} | {25.0:<10.2f} | {guess_comp.b:<10.2f} | {fitted_comp.b:<10.2f}")
    print("-" * 30)
    
    return spectrum, fitter, new_syst_list, result


def _plot_fit(spectrum, fitter, new_syst_list, result):
    """Plot the fit result (only when run as a script)."""
    try:
        x_plot, y_model = fitter.compute_model_flux()
        # We need to re-compute model for the fitted result manually to be sure
        # (Fitter.compute_model_flux uses the *current* constraints, but the fitter returns a *new* list)
        
        # To plot the result, we make a new fitter with the new list (statelessness check)
        fitter_final = VoigtFitter(spectrum, new_syst_list)
        x_plot, y_final = fitter_final.compute_model_flux()
        
        plt.figure(figsize=(10, 6))
        plt.step(spectrum.x.value, spectrum.y.value, color='black', alpha=0.3, label='Data (Noisy)')
        plt.plot(x_plot, y_final, color='red', lw=2, label='Fitted Model')
        plt.title("V2 Voigt Fitter Test: CIV Doublet")
        plt.xlabel("Wavelength (Angstrom)")
        plt.ylabel("Normalized Flux")
        plt.legend()
        plt.show()
        print("Plot displayed.")
    except Exception as e:
        print(f"Plotting skipped: {e}")

def test_fitting_workflow():
    _, _, new_syst_list, result = _run_fit()
    fitted = new_syst_list.components[0]

    # The fit starts from a wrong guess and must recover the true parameters
    # (z=2.0, logN=14.0, b=25.0) within the noise.
    assert result.success
    assert abs(fitted.z - 2.0) < 1e-4
    assert abs(fitted.logN - 14.0) < 0.1
    assert abs(fitted.b - 25.0) < 3.0


if __name__ == "__main__":
    _plot_fit(*_run_fit())