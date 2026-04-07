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

# Keep logs quiet so we can see the report
logging.basicConfig(level=logging.WARNING, format='[%(levelname)s] %(message)s')

snr = 200.0

def make_stress_spectrum():
    print("--- GENERATING 'TRAIN WRECK' SYNTHETIC DATA ---")
    x = np.linspace(154.5, 155.5, 2000) 
    cont = np.ones_like(x)
    
    # 10 Components (Ground Truth)
    truth = [
        # Main Clump
        {'z': 0.00000, 'logN': 14.0, 'b': 12.0},
        {'z': 0.00005, 'logN': 13.8, 'b': 10.0},
        {'z': -0.00005,'logN': 13.5, 'b': 8.0},  
        # Blue Wing
        {'z': -0.00040, 'logN': 13.2, 'b': 15.0},
        {'z': -0.00050, 'logN': 13.5, 'b': 10.0},
        {'z': -0.00060, 'logN': 13.0, 'b': 6.0}, 
        # High Velocity
        {'z': 0.00100, 'logN': 13.7, 'b': 20.0},
        {'z': 0.00110, 'logN': 13.3, 'b': 10.0},
        # Background
        {'z': 0.00020, 'logN': 12.8, 'b': 25.0}, 
        {'z': -0.00020,'logN': 12.9, 'b': 22.0}
    ]
    for t in truth: t['series'] = 'CIV'

    systs_data = SystemListDataV2()
    gt_systs = SystemListV2(systs_data)
    for c in truth:
        gt_systs = gt_systs.add_component(**c)

    # Note: R=50000 -> FWHM ~ 6 km/s
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
    
    rng = np.random.default_rng(999) 
    noise = rng.normal(0, 1/snr, len(x)) 
    flux_noisy = flux_true + noise
    dy = np.ones_like(x) / snr
    
    real_spec_data = SpectrumDataV2(
        x=DataColumnV2(x, au.nm),
        xmin=DataColumnV2(x, au.nm), xmax=DataColumnV2(x, au.nm),
        y=DataColumnV2(flux_noisy, au.dimensionless_unscaled),
        dy=DataColumnV2(dy, au.dimensionless_unscaled),
        aux_cols={'cont': DataColumnV2(cont, au.dimensionless_unscaled)},
        z_em=0.1, z_rf=0.0, resol=50000.0, meta={}
    )
    return SpectrumV2(real_spec_data), truth

def run_stress_test():
    spec, truth_list = make_stress_spectrum()
    
    print("\n--- INITIALIZING FIT ---")
    systs = SystemListV2(SystemListDataV2())
    systs = systs.add_component(series='CIV', z=0.0, logN=14.0, b=20.0) 
    seed_uuid = systs.components[0].uuid
    
    print("--- RUNNING STRESS OPTIMIZATION ---")
    # Using min_dv=5.0 to balance de-blending vs stability
    final_systs = systs.optimize_hierarchy(
        spec=spec,
        uuid_seed=seed_uuid,
        max_components=15, 
        threshold_sigma=1.0, 
        aic_penalty=1.0,     
        z_window_kms=600.0,
        min_dv=5.0,          # <-- TUNED: 2.0 might be too unstable, 10.0 too coarse
        patience=3           
    )
    
    # --- REPORTING ---
    print("\n" + "="*40)
    print("          VALIDATION REPORT")
    print("="*40)
    
    comps = final_systs.components
    print(f"Ground Truth: {len(truth_list)} components")
    print(f"Recovered:    {len(comps)} components")
    
    # Matching Logic
    unmatched_truth = truth_list.copy()
    matches = []
    
    print("\n[MATCHED COMPONENTS]")
    print(f"{'d_vel':<8} | {'z_fit':<9} | {'logN (Fit/True)':<18} | {'b (Fit/True)':<15}")
    print("-" * 65)
    
    # Sort by position for cleaner output
    sorted_comps = sorted(comps, key=lambda c: c.z)
    
    c_kms = 299792.458
    
    for c_rec in sorted_comps:
        best_dist = 1e9
        best_t_idx = -1
        
        # Find closest truth
        for i, t in enumerate(unmatched_truth):
            dist = abs(c_rec.z - t['z'])
            if dist < best_dist:
                best_dist = dist
                best_t_idx = i
        
        # Match threshold: 30 km/s
        if best_dist < 0.0003: 
            t = unmatched_truth.pop(best_t_idx)
            dv = c_kms * (c_rec.z - t['z']) / (1 + t['z'])
            print(f"{dv:6.1f} k | {c_rec.z:9.5f} | {c_rec.logN:5.2f} / {t['logN']:5.2f}   | {c_rec.b:4.1f} / {t['b']:4.1f}")
        else:
            print(f"{'---':<8} | {c_rec.z:9.5f} | {c_rec.logN:5.2f} / {'---':<5}   | {c_rec.b:4.1f} / {'---':<4} (Spurious?)")

    print("\n[MISSED TRUTHS]")
    if not unmatched_truth:
        print("None! Perfect Recall.")
    for t in unmatched_truth:
        print(f"Missed z={t['z']:.5f} (logN={t['logN']}, b={t['b']})")

    print("="*40)

    # --- PLOTTING ---
    print("\nLaunching Plotter...")
    fitter = VoigtFitterV2(spec, final_systs)
    _, model_flux = fitter.compute_model_flux()
    
    x = spec.x.to(au.nm).value
    y = spec.y.value
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
    ax1.step(x, y, where='mid', color='grey', alpha=0.4, label='Noisy Data')
    ax1.plot(x, model_flux, color='blue', lw=1.5, label='Fit Model')
    ax1.plot(x, spec.cont.value, color='orange', ls=':', alpha=0.5)
    
    lam_0 = ATOM_DATA['CIV_1548']['wave']/10.0
    
    # Draw Truth (Green)
    for t in truth_list:
        lx = lam_0*(1+t['z'])
        ax1.axvline(lx, color='green', alpha=0.3, ymax=0.15)
        
    # Draw Recovered (Red)
    for c in sorted_comps:
        lx = lam_0*(1+c.z)
        ax1.axvline(lx, color='red', alpha=0.8, ymin=0.0, ymax=0.2)
        ax1.text(lx, 0.2, f"{c.b:.0f}", color='red', fontsize=8, rotation=90, transform=ax1.get_xaxis_transform())

    ax1.set_title(f"Stress Test (SNR={snr}, min_dv=5.0)\nRed numbers = Fitted b-parameter")
    ax1.legend()
    
    resid = (y - model_flux) / spec.dy.value
    ax2.step(x, resid, where='mid', color='black', lw=0.5)
    ax2.axhline(2, color='red', ls=':', alpha=0.3)
    ax2.axhline(-2, color='red', ls=':', alpha=0.3)
    ax2.set_ylim(-5, 5)
    ax2.set_ylabel("Residuals")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_stress_test()