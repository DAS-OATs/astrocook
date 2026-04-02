import numpy as np
import astropy.units as au
import matplotlib.pyplot as plt
import logging
import argparse
import sys
import os

# Import Astrocook Core
from astrocook.core.structures import (
    SpectrumDataV2, DataColumnV2, ComponentDataV2, SystemListDataV2
)
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.system_list import SystemListV2
from astrocook.fitting.voigt_fitter import VoigtFitterV2
from astrocook.fitting.bayesian_fitter import BayesianVoigtFitter

def generate_synthetic_data(n_components=1, snr=20, regime='unsaturated'):
    """
    Generates a synthetic spectrum with tunable complexity.
    
    Regimes:
    - unsaturated: logN ~ 13-14 (linear region)
    - saturated: logN ~ 15-16 (flat region of CoG)
    - dla: logN ~ 20 (damping wings)
    """
    logging.info(f"Generating data: {n_components} components, SNR={snr}, regime={regime}")
    
    z_center = 2.0
    lambda_0 = 1215.67  # Ly-alpha
    lambda_obs = lambda_0 * (1 + z_center)
    
    # Grid
    x_arr = np.linspace(lambda_obs - 20, lambda_obs + 20, 1000)
    
    components = []
    for i in range(n_components):
        # Slightly offset components
        z_offset = (i - (n_components-1)/2.0) * 0.0005
        z_i = z_center + z_offset
        
        if regime == 'unsaturated':
            logN_i = 13.5 + np.random.uniform(-0.2, 0.2)
            b_i = 20.0 + np.random.uniform(-5, 5)
        elif regime == 'saturated':
            logN_i = 15.5 + np.random.uniform(-0.5, 0.5)
            b_i = 15.0 + np.random.uniform(-5, 5)
        elif regime == 'dla':
            logN_i = 20.3
            b_i = 30.0
        
        comp = ComponentDataV2(
            id=i, z=z_i, dz=None, logN=logN_i, dlogN=None, b=b_i, db=None, series='Ly_a'
        )
        components.append(comp)
        logging.info(f"  Component {i}: z={z_i:.5f}, logN={logN_i:.2f}, b={b_i:.1f}")

    sys_list = SystemListV2(SystemListDataV2(components=components))
    
    # Create dummy spectrum for generation
    x_col = DataColumnV2(x_arr, au.Angstrom)
    spec_data = SpectrumDataV2(
        x=x_col, xmin=x_col, xmax=x_col, 
        y=DataColumnV2(np.ones_like(x_arr), au.dimensionless_unscaled),
        dy=DataColumnV2(np.ones_like(x_arr)/snr, au.dimensionless_unscaled)
    )
    spec = SpectrumV2(spec_data)
    
    # Generate model
    fitter = VoigtFitterV2(spec, sys_list)
    # We use Case B loop for generation by not setting up fit context
    _, true_flux_phys = fitter.compute_model_flux()
    true_flux = true_flux_phys / fitter._norm
    
    # Add Noise
    noise = np.random.normal(0, 1.0/snr, size=len(x_arr))
    y_obs = true_flux + noise
    dy_obs = np.ones_like(y_obs) / snr
    
    # Final Spectrum
    final_data = SpectrumDataV2(
        x=x_col, xmin=x_col, xmax=x_col,
        y=DataColumnV2(y_obs, au.dimensionless_unscaled),
        dy=DataColumnV2(dy_obs, au.dimensionless_unscaled),
        aux_cols={'cont': DataColumnV2(np.ones_like(x_arr), au.dimensionless_unscaled)}
    )
    return SpectrumV2(final_data), sys_list

def run_comparison(n_components=1, snr=20, regime='unsaturated'):
    # 1. Setup
    spec, true_sys = generate_synthetic_data(n_components, snr, regime)
    
    # 2. Frequentist Fit
    logging.info("\n--- Running Frequentist Fit ---")
    # Start with a slightly perturbed guess
    guess_comps = []
    for c in true_sys.components:
        gc = ComponentDataV2(
            id=c.id, z=c.z + 0.0001, dz=None, logN=c.logN - 0.5, dlogN=None, b=c.b + 5.0, db=None, series=c.series
        )
        guess_comps.append(gc)
    guess_sys = SystemListV2(SystemListDataV2(components=guess_comps))
    
    f_fitter = VoigtFitterV2(spec, guess_sys)
    fit_sys, model_flux, res = f_fitter.fit(verbose=0)
    
    logging.info(f"Frequentist Success: {res.success}")
    for i, (t, f) in enumerate(zip(true_sys.components, fit_sys.components)):
        logging.info(f"Comp {i} - logN: True={t.logN:.2f}, Fit={f.logN:.2f} +/- {f.dlogN:.2f}")
    
    # 3. Bayesian Fit
    logging.info("\n--- Running Bayesian Fit (MCMC) ---")
    # NEW API: BayesianVoigtFitter(spec, systs)
    b_fitter = BayesianVoigtFitter(spec, guess_sys)
    
    try:
        mcmc_res = b_fitter.run_mcmc(nwalkers=32, nsteps=500, burn_in=100)
        medians = mcmc_res['medians']
        stds = mcmc_res['std']
        
        logging.info("Bayesian (MCMC) Results:")
        logging.info(f"  MCMC Medians: {medians}")
        logging.info(f"  MCMC StdDevs: {stds}")
        
    except Exception as e:
        logging.error(f"MCMC fit failed: {e}")
        mcmc_res = None

    # 4. Bayesian Fit (Nested Sampling)
    logging.info("\n--- Running Bayesian Fit (Nested Sampling) ---")
    try:
        nested_res = b_fitter.run_nested(nlive=250)
        logging.info(f"Nested logZ: {nested_res['logz']:.2f} +/- {nested_res['logzerr']:.2f}")
    except Exception as e:
        logging.error(f"Nested Sampling failed: {e}")
        nested_res = None

    # 5. Visualization (Corner Plot)
    if mcmc_res is not None:
        try:
            import corner
            # Use samples from MCMC for corner plot
            labels = ["z", "logN", "b"]
            fig = corner.corner(mcmc_res['samples'], labels=labels, truths=[true_sys.components[0].z, true_sys.components[0].logN, true_sys.components[0].b])
            plot_path = os.path.abspath("tests/bayesian_corner_plot.png")
            fig.savefig(plot_path)
            logging.info(f"Corner plot saved to {plot_path}")
        except ImportError:
            logging.warning("corner library not found, skipping plot.")
        except Exception as e:
            logging.error(f"Corner plot failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--components", type=int, default=1)
    parser.add_argument("--snr", type=float, default=20)
    parser.add_argument("--regime", type=str, default='unsaturated')
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    run_comparison(args.components, args.snr, args.regime)
