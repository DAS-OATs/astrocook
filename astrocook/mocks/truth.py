"""Ground truth of a mock spectrum."""
import dataclasses
import json
from typing import List, Optional


@dataclasses.dataclass(frozen=True)
class MockComponent:
    """One Voigt component.

    Parameters
    ----------
    series : str
        A transition (e.g. ``'CIV_1548'``), a multiplet of
        ``STANDARD_MULTIPLETS`` (e.g. ``'CIV'``) or a comma-separated list
        of transitions.
    z : float
        Redshift.
    logN : float
        Log10 of the column density (cm^-2).
    b : float
        Doppler parameter (km/s).
    btur : float, optional
        Turbulent broadening (km/s), added in quadrature to ``b``.
    """
    series: str
    z: float
    logN: float
    b: float
    btur: float = 0.0


@dataclasses.dataclass(frozen=True)
class MockTruth:
    """Everything needed to generate, and later evaluate, a mock spectrum.

    The wavelength grid is uniform in velocity (pixel size ``dv_kms``), as
    for echelle spectrographs. The continuum is a power law
    ``cont_amplitude * (lambda / cont_ref_nm) ** cont_slope``. The noise is
    Gaussian, with a constant signal-to-noise ratio ``snr`` per pixel at the
    continuum level (``None``: noiseless).

    ``resol`` is the resolving power, either a single value or ``None`` (no
    instrumental broadening). ``oversample`` is the number of sub-pixels on
    which the model is computed before being averaged over each pixel.
    """
    components: List[MockComponent]
    lam_min_nm: float
    lam_max_nm: float
    dv_kms: float
    resol: Optional[float] = 45000.0
    snr: Optional[float] = 50.0
    z_em: float = 0.0
    seed: int = 0
    oversample: int = 5
    cont_amplitude: float = 1.0
    cont_slope: float = 0.0
    cont_ref_nm: float = 500.0

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "MockTruth":
        data = dict(data)
        data["components"] = [MockComponent(**c) for c in data["components"]]
        return cls(**data)

    def to_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_json(cls, path: str) -> "MockTruth":
        with open(path) as f:
            return cls.from_dict(json.load(f))
