import numpy as np
import astropy.units as au
import logging
import matplotlib.pyplot as plt

# Import V2 Core
from astrocook.core.structures import (
    SpectrumDataV2, DataColumnV2, ComponentDataV2, SystemListDataV2
)
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.core.atomic_data import ATOM_DATA

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
    true_comp = ComponentDataV2(
        id=0, z=z_true, dz=0, logN=14.0, dlogN=0, b=25.0, db=0, 
        series='CIV' # This triggers both 1548 and 1550
    )
    
    # Create a dummy system list for generation
    dummy_list_data = SystemListDataV2(components=[true_comp])
    dummy_sys = SystemListV2(dummy_list_data)
    
    # Create a dummy spectrum container
    x_col = DataColumnV2(x_arr, au.Angstrom)
    y_dummy = DataColumnV2(np.ones_like(x_arr), au.dimensionless_unscaled)
    dy_dummy = DataColumnV2(np.ones_like(x_arr)*0.1, au.dimensionless_unscaled)
    spec_data = SpectrumDataV2(x=x_col, xmin=x_col, xmax=x_col, y=y_dummy, dy=dy_dummy)
    spec_v2 = SpectrumV2(spec_data)
    
    # Use Fitter to generate true flux
    fitter = VoigtFitterV2(spec_v2, dummy_sys)
    # We need to hack the constraint model to give us the P vector for the true component
    # Since we are just generating, we can calculate manually or use the internal method if accessible
    # Let's trust the fitter logic:
    # We create a fitter, but we need to generate the model from the *parameters*, not fit.
    # VoigtFitterV2._compute_model uses the constraints vector.
    
    # Let's extract the p_vector corresponding to the true component
    p_true = np.array([z_true, 14.0, 25.0, 0.0]) # z, logN, b, btur
    true_flux = fitter._compute_model(p_true)
    
    # 3. Add Noise
    noise_level = 0.05
    noise = np.random.normal(0, noise_level, size=len(x_arr))
    y_obs = true_flux + noise
    dy_obs = np.ones_like(y_obs) * noise_level
    
    # 4. Return SpectrumV2 Object
    y_col = DataColumnV2(y_obs, au.dimensionless_unscaled)
    dy_col = DataColumnV2(dy_obs, au.dimensionless_unscaled)
    
    # Add a continuum column (all ones)
    aux = {'cont': DataColumnV2(np.ones_like(x_arr), au.dimensionless_unscaled)}
    
    final_data = SpectrumDataV2(
        x=x_col, xmin=x_col, xmax=x_col, y=y_col, dy=dy_col, aux_cols=aux
    )
    return SpectrumV2(final_data)

def test_fitting_workflow():
    
    # 1. Get Data
    spectrum = create_synthetic_spectrum()
    
    # 2. Create a "Bad" Initial Guess (User clicked slightly wrong)
    # True was z=2.0, logN=14.0, b=25.0
    # Guess is z=2.001, logN=13.5, b=40.0
    print("\n--- Initializing System List with Guess ---")
    guess_comp = ComponentDataV2(
        id=1, 
        z=2.0005,  # Slightly shifted
        dz=None,
        logN=13.5, # Too weak
        dlogN=None,
        b=40.0,    # Too broad
        db=None,
        series='CIV'
    )
    
    syst_list = SystemListV2(SystemListDataV2(components=[guess_comp]))
    
    print(f"Guess: z={guess_comp.z:.5f}, logN={guess_comp.logN:.2f}, b={guess_comp.b:.2f}")
    
    # 3. Initialize Fitter
    print("\n--- Initializing Fitter ---")
    fitter = VoigtFitterV2(spectrum, syst_list)
    
    # 4. Run Fit
    print("--- Running Fit (Levenberg-Marquardt) ---")
    new_syst_list, result = fitter.fit(verbose=2)
    
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
    
    # 6. Visualization (Optional)
    try:
        x_plot, y_model = fitter.compute_model_flux()
        # We need to re-compute model for the fitted result manually to be sure
        # (Fitter.compute_model_flux uses the *current* constraints, but the fitter returns a *new* list)
        
        # To plot the result, we make a new fitter with the new list (statelessness check)
        fitter_final = VoigtFitterV2(spectrum, new_syst_list)
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

if __name__ == "__main__":
    test_fitting_workflow()