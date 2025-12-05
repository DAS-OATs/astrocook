import sys
import os
import logging
import numpy as np
import astropy.units as au
import matplotlib.pyplot as plt

# Ensure import path is correct
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from astrocook.core.structures import SpectrumDataV2, DataColumnV2, SystemListDataV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.core.atomic_data import ATOM_DATA

# Configure Logging
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

def make_synthetic_spectrum():
    print("\n--- GENERATING SYNTHETIC DATA ---")
    x = np.linspace(154.8, 155.5, 1000) # nm (CIV region)
    cont = np.ones_like(x)
    
    truth = [
        {'series': 'CIV', 'z': 0.0,      'logN': 14.0, 'b': 15.0},
        {'series': 'CIV', 'z': 0.0005,   'logN': 13.0, 'b': 10.0}
    ]
    
    systs_data = SystemListDataV2()
    gt_systs = SystemListV2(systs_data)
    for c in truth:
        gt_systs = gt_systs.add_component(**c)
        print(f"Truth: {c['series']} at z={c['z']:.5f} (logN={c['logN']}, b={c['b']})")

    mock_spec_data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x, au.nm), xmax=DataColumnV2(x, au.nm),
        y=DataColumnV2(cont, au.dimensionless_unscaled),
        dy=DataColumnV2(np.ones_like(x)*0.1, au.dimensionless_unscaled),
        aux_cols={'cont': DataColumnV2(cont, au.dimensionless_unscaled)},
        z_em=0.1, z_rf=0.0, resol=50000.0, meta={}
    )
    mock_spec = SpectrumV2(mock_spec_data)
    
    fitter = VoigtFitterV2(mock_spec, gt_systs)
    _, flux_true = fitter.compute_model_flux()
    
    rng = np.random.default_rng(42) 
    noise = rng.normal(0, 0.05, len(x))
    flux_noisy = flux_true + noise
    dy = np.ones_like(x) * 0.05
    
    real_spec_data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x, au.nm), xmax=DataColumnV2(x, au.nm),
        y=DataColumnV2(flux_noisy, au.dimensionless_unscaled),
        dy=DataColumnV2(dy, au.dimensionless_unscaled),
        aux_cols={'cont': DataColumnV2(cont, au.dimensionless_unscaled)},
        z_em=0.1, z_rf=0.0, resol=50000.0, meta={}
    )
    return SpectrumV2(real_spec_data), truth

def plot_results(spec, final_systs, truth_list):
    """Visualizes the fit against truth."""
    print("\n--- PLOTTING ---")
    
    # 1. Compute Model
    fitter = VoigtFitterV2(spec, final_systs)
    _, model_flux = fitter.compute_model_flux()
    
    x = spec.x.to(au.nm).value
    y = spec.y.value
    dy = spec.dy.value
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
    # PANEL 1: Flux
    ax1.step(x, y, where='mid', color='grey', alpha=0.5, label='Data')
    ax1.plot(x, model_flux, color='blue', linewidth=1.5, label='Fit Model')
    
    # Plot Truth
    for t in truth_list:
        lam_0 = ATOM_DATA['CIV_1548']['wave']/10.0
        lam_obs = lam_0 * (1 + t['z'])
        ax1.axvline(lam_obs, color='green', linestyle='-', alpha=0.6, ymax=0.2, label='Truth' if t==truth_list[0] else "")
        
    # Plot Recovered
    for c in final_systs.components:
        if '1548' in c.series or c.series == 'CIV':
            lam_0 = ATOM_DATA['CIV_1548']['wave']/10.0
            lam_obs = lam_0 * (1 + c.z)
            ax1.axvline(lam_obs, color='red', linestyle='--', alpha=0.8, ymax=0.2, label='Recovered' if c==final_systs.components[0] else "")

    ax1.set_ylabel("Normalized Flux")
    ax1.legend()
    ax1.set_title("Standard Verification Test (2 Components)")

    # PANEL 2: Residuals
    resid = (y - model_flux) / dy
    ax2.step(x, resid, where='mid', color='black', lw=0.8)
    ax2.axhline(0, color='red', alpha=0.5)
    ax2.axhline(2, color='red', linestyle=':', alpha=0.3)
    ax2.axhline(-2, color='red', linestyle=':', alpha=0.3)
    ax2.set_ylabel("Residuals (sigma)")
    ax2.set_xlabel("Wavelength (nm)")
    ax2.set_ylim(-5, 5)
    
    plt.tight_layout()
    plt.show()

def run_test():
    spec, truth_list = make_synthetic_spectrum()
    
    print("\n--- INITIALIZING FIT ---")
    systs = SystemListV2(SystemListDataV2())
    systs = systs.add_component(series='CIV', z=0.0, logN=13.8, b=20.0) 
    seed_uuid = systs.components[0].uuid
    
    print("\n--- RUNNING OPTIMIZE_HIERARCHY ---")
    final_systs = systs.optimize_hierarchy(
        spec=spec,
        uuid_seed=seed_uuid,
        max_components=3,
        threshold_sigma=2.0,
        aic_penalty=0.0, 
        z_window_kms=300.0,
        patience=2
    )
    
    # ... (Keep existing console print logic) ...
    print("\n--- RESULTS ---")
    comps = final_systs.components
    print(f"Recovered {len(comps)} components (Expected {len(truth_list)})")
    
    # [VISUALIZATION]
    plot_results(spec, final_systs, truth_list)

if __name__ == "__main__":
    run_test()