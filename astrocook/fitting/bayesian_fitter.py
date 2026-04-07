import numpy as np
import logging
from typing import Tuple, Dict, Any, Optional
from .voigt_fitter import VoigtFitterV2

try:
    import emcee
except ImportError:
    emcee = None

try:
    import dynesty
except ImportError:
    dynesty = None

class BayesianVoigtFitter:
    """
    Bayesian engine for Voigt profile fitting.
    
    This class wraps a :class:`~astrocook.fitting.voigt_fitter.VoigtFitterV2`
    instance to perform Markov Chain Monte Carlo (MCMC) or Nested Sampling
    using the exact same physical model, constraints, and data preparation.

    Parameters
    ----------
    spec : SpectrumV2
        The spectral data to fit.
    systs : SystemListV2
        The system list containing the components to model.
    z_window_kms : float, optional
        Velocity window (km/s) around component centers to include in the fit mask.
        Defaults to ``20.0``.
    """
    def __init__(self, spec: Any, systs: Any, z_window_kms: float = 20.0):
        self.spec = spec
        self.systs = systs
        self.z_window_kms = z_window_kms
        
        # We need a VoigtFitterV2 instance to handle the model and constraints
        self.fitter = VoigtFitterV2(self.spec, self.systs)
        
        # Prepare the data context from the standard fitter
        self.lower_bounds, self.upper_bounds = self.fitter.prepare_fit_context(z_window_kms=z_window_kms)
        self.ctx = self.fitter.get_fit_context_data()
        
        self.p_initial = self.fitter._constraints.p_free_vector
        self.ndim = len(self.p_initial)
        
    def ln_likelihood(self, p_free: np.ndarray) -> float:
        """
        Compute the log-likelihood of a parameter vector.

        Evaluates the standard Gaussian log-likelihood for the active fit region,
        incorporating any penalty terms (e.g., Doppler parameter runaways) 
        inherited from the frequentist fitter.

        Parameters
        ----------
        p_free : np.ndarray
            The vector of free parameters to evaluate.

        Returns
        -------
        float
            The computed log-likelihood, or ``-np.inf`` if evaluation fails.
        """
        try:
            model = self.ctx['compute_model'](p_free)
            # Only use pixels in the dynamic mask
            mask = self.ctx['mask']
            resid = (self.ctx['y'][mask] - model[mask]) / self.ctx['dy'][mask]
            
            # Additional penalty terms from Frequentist fitter (e.g. b-parameter runaway)
            # To match the logic of VoigtFitterV2._residual_function
            p_full = self.fitter._constraints.map_p_free_to_full(p_free)
            is_free_mask = self.fitter._constraints._param_map['is_free']
            num_params = 4
            penalty_chi2 = 0.0
            
            for i in range(len(p_full) // num_params):
                idx_b = i * num_params + 2
                if is_free_mask[idx_b]:
                    b_val = p_full[idx_b]
                    if b_val > 50.0:
                        penalty_chi2 += (5.0 * (b_val - 50.0))**2
            
            return -0.5 * (np.sum(resid**2) + penalty_chi2)
        except Exception as e:
            return -np.inf

    def ln_prior(self, p_free: np.ndarray) -> float:
        """
        Compute the log-prior probability of a parameter vector.

        Applies a uniform flat prior defined by the physical boundaries of the
        free parameters (e.g., column density between 10.0 and 23.0).

        Parameters
        ----------
        p_free : np.ndarray
            The vector of free parameters.

        Returns
        -------
        float
            ``0.0`` if the parameters are within physical bounds, ``-np.inf`` otherwise.
        """
        if np.any(p_free < self.lower_bounds) or np.any(p_free > self.upper_bounds):
            return -np.inf
        return 0.0 # Log of 1/Volume is constant, so we can use 0 for simplicity

    def ln_posterior(self, p_free: np.ndarray) -> float:
        """
        Compute the total log-posterior probability.

        Calculates the sum of the log-prior and log-likelihood.

        Parameters
        ----------
        p_free : np.ndarray
            The vector of free parameters.

        Returns
        -------
        float
            The log-posterior probability, or ``-np.inf`` if outside prior bounds.
        """
        lp = self.ln_prior(p_free)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.ln_likelihood(p_free)

    def run_mcmc(self, nwalkers: int = 32, nsteps: int = 1000, 
                 burn_in: int = 200, progress: bool = True) -> Dict[str, Any]:
        """
        Execute Markov Chain Monte Carlo (MCMC) sampling using ``emcee``.

        Initializes walkers around the initial parameter guess and samples the
        posterior distribution.

        Parameters
        ----------
        nwalkers : int, optional
            Number of MCMC walkers. Defaults to ``32``.
        nsteps : int, optional
            Number of sampling steps per walker. Defaults to ``1000``.
        burn_in : int, optional
            Number of initial steps to discard as burn-in. Defaults to ``200``.
        progress : bool, optional
            If True, pipes progress percentage to the logger. Defaults to ``True``.

        Returns
        -------
        dict
            A dictionary containing:
            - ``'samples'``: The flattened, burn-in discarded chain (np.ndarray).
            - ``'sampler'``: The raw ``emcee.EnsembleSampler`` object.
            - ``'medians'``: The median values for each parameter (np.ndarray).
            - ``'std'``: The standard deviations for each parameter (np.ndarray).

        Raises
        ------
        ImportError
            If the ``emcee`` package is not installed.
        """
        if emcee is None:
            raise ImportError("emcee is not installed. Please install it to use MCMC.")
            
        # Initialize walkers in a small cluster around the best-fit / initial guess
        # We add some noise but stay within bounds
        pos = []
        for _ in range(nwalkers):
            p = self.p_initial + 1e-4 * np.random.randn(self.ndim)
            p = np.clip(p, self.lower_bounds + 1e-8, self.upper_bounds - 1e-8)
            pos.append(p)
        pos = np.array(pos)
        
        sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.ln_posterior)
        
        logging.info(f"Running MCMC: {nwalkers} walkers, {nsteps} steps...")
        if progress:
            # We use manual iteration instead of run_mcmc to pipe progress to the GUI via logging
            # Astrocook's RecipeProgressDialog looks for '>> PROGRESS: X' messages
            update_interval = max(1, nsteps // 100)

            for i, result in enumerate(sampler.sample(pos, iterations=nsteps, progress=False)):
                
                # Update GUI every 1% or at least every step
                if i % update_interval == 0:
                    pct = int(100 * i / nsteps)
                    logging.info(f">> PROGRESS: {pct}")
            
            logging.info(">> PROGRESS: 100")
        else:
            sampler.run_mcmc(pos, nsteps, progress=False)
        
        # Process results
        flat_samples = sampler.get_chain(discard=burn_in, flat=True)
        
        return {
            'samples': flat_samples,
            'sampler': sampler,
            'medians': np.median(flat_samples, axis=0),
            'std': np.std(flat_samples, axis=0)
        }

    def run_nested(self, nlive: int = 250, progress: bool = True) -> Dict[str, Any]:
        """
        Execute Nested Sampling using ``dynesty``.

        Samples the parameter space to estimate both the posterior distribution
        and the Bayesian evidence (log-Z).

        Parameters
        ----------
        nlive : int, optional
            Number of live points for the nested sampler. Defaults to ``250``.
        progress : bool, optional
            If True, pipes progress (dLogZ) to the logger. Defaults to ``True``.

        Returns
        -------
        dict
            A dictionary containing:
            - ``'results'``: The raw ``dynesty`` results object.
            - ``'logz'``: The estimated log-evidence.
            - ``'logzerr'``: The error on the log-evidence.
            - ``'samples'``: The posterior samples (np.ndarray).
            - ``'medians'``: The median values for each parameter (np.ndarray).
            - ``'std'``: The standard deviations for each parameter (np.ndarray).

        Raises
        ------
        ImportError
            If the ``dynesty`` package is not installed.
        """
        if dynesty is None:
            raise ImportError("dynesty is not installed. Please install it to use Nested Sampling.")
            
        def prior_transform(u):
            """Transform uniform unit cube to physical bounds."""
            return self.lower_bounds + (self.upper_bounds - self.lower_bounds) * u

        dsampler = dynesty.NestedSampler(self.ln_likelihood, prior_transform, self.ndim, nlive=nlive)
        
        logging.info(f"Running Nested Sampling (nlive={nlive})...")
        if progress:
            import time
            last_log_time = time.time()
            
            # Start manual sampling to safely pipe progress to GUI
            # We target dlogz=0.1 which is dynesty's standard default
            target_dlogz = 0.1
            for i, res in enumerate(dsampler.sample(dlogz=target_dlogz)):
                if time.time() - last_log_time > 0.5:
                    dlogz = getattr(dsampler, 'dlogz', None)
                    if dlogz is None and isinstance(res, tuple) and len(res) > 3:
                        dlogz = res[-1]  # fallback for some old dynesty versions
                    
                    if dlogz is not None:
                        logging.info(f"Nested Iter {i} | dLogZ: {dlogz:.2f} > {target_dlogz}")
                    else:
                        logging.info(f"Nested Iter {i}...")
                    last_log_time = time.time()
            
            # Finalize remaining live points
            logging.info("Adding final live points...")
            dsampler.add_live_points()
        else:
            dsampler.run_nested(print_progress=False)
        
        res = dsampler.results
        samples = res.samples
        
        return {
            'results': res,
            'logz': res.logz[-1],
            'logzerr': res.logzerr[-1],
            'samples': samples,
            'medians': np.median(samples, axis=0),
            'std': np.std(samples, axis=0)
        }

    def apply_results(self, res: Dict[str, Any]) -> Any:
        """
        Apply Bayesian medians and uncertainties back to the system list.

        Extracts the posterior medians and standard deviations, maps them back
        to the full parameter vector, and updates the corresponding components
        in the system list. Fixed or linked parameters maintain their prior state.

        Parameters
        ----------
        res : dict
            The results dictionary returned by :meth:`run_mcmc` or :meth:`run_nested`.

        Returns
        -------
        SystemListV2
            A new system list containing the updated parameters and Bayesian errors.
        """
        medians = res['medians'] if 'medians' in res else np.median(res['samples'], axis=0)
        stds = res['std'] if 'std' in res else np.std(res['samples'], axis=0)
        
        # Map p_free back to full parameter vector
        p_full = self.fitter._constraints.map_p_free_to_full(medians)
        
        # Map uncertainties back to full vector (treating frozen as 0 uncertainty)
        dp_full = np.zeros_like(p_full)
        is_free_mask = self.fitter._constraints._param_map['is_free']
        dp_full[is_free_mask] = stds
        
        # We only update components that were part of the active group
        # This ensures we don't accidentally reset/clear errors for other components
        new_systs = self.systs
        active_uuids = self.fitter._constraints._active_uuids
        if active_uuids is None:
            # If everything was active (no context), update everything
            active_uuids = [c.uuid for c in self.fitter._constraints._system_list.components]

        for uuid in active_uuids:
            # Find its offset in p_full
            idx = self.fitter._constraints._uuid_to_idx.get(uuid)
            if idx is not None:
                # Update basic parameters (medians)
                updates = {
                    'z': p_full[idx],
                    'logN': p_full[idx+1],
                    'b': p_full[idx+2],
                    'btur': p_full[idx+3]
                }
                
                # Update uncertainties ONLY for parameters that were free
                # Non-free parameters (Frozen or Linked) should keep their previous errors
                # or stay at None if they didn't have any.
                if is_free_mask[idx]: updates['dz'] = dp_full[idx]
                if is_free_mask[idx+1]: updates['dlogN'] = dp_full[idx+1]
                if is_free_mask[idx+2]: updates['db'] = dp_full[idx+2]
                if is_free_mask[idx+3]: updates['dbtur'] = dp_full[idx+3]
                
                new_systs = new_systs.update_component(uuid, **updates)
        
        return new_systs
