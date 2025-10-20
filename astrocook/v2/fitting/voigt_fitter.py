class VoigtFitterV2:
    """Orchestrates the minimization using SciPy, PyMC, or other backends."""
    
    def __init__(self, spectrum: SpectrumV2, constraint_model: ConstraintModelV2):
        self.spectrum = spectrum
        self.constraints = constraint_model

    def fit(self, method: str = 'least_squares'):
        # 1. Get P_free vector
        p_free = self.constraints.p_free_vector
        
        # 2. Select minimization backend
        if method == 'least_squares':
            # Calls SciPy.optimize.least_squares using a pure residual function
            pass
        elif method == 'bayesian':
            # Calls emcee or PyMC sampler
            pass
        
        # 3. Return the fit result and a new immutable SystemListV2
        pass # Implementation needed