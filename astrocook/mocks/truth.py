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
class EmissionLine:
    """Gaussian emission line added to the continuum (quasar frame).

    Parameters
    ----------
    rest_nm : float
        Rest-frame wavelength (nm); observed at ``rest_nm * (1 + z_em)``.
    amplitude : float
        Peak height, in units of the power-law continuum.
    fwhm_kms : float
        Full width at half maximum (km/s).
    """
    rest_nm: float
    amplitude: float
    fwhm_kms: float


#: Lyman-alpha + NV and CIV emission, indicative values
QSO_EMISSION_LINES = [EmissionLine(121.567, 2.0, 4000.0),
                      EmissionLine(123.9, 0.4, 3000.0),
                      EmissionLine(154.9, 0.8, 5000.0)]


@dataclasses.dataclass(frozen=True)
class MockTruth:
    """Everything needed to generate, and later evaluate, a mock spectrum.

    The wavelength grid is uniform in velocity (pixel size ``dv_kms``), as
    for echelle spectrographs. The continuum is a power law
    ``cont_amplitude * (lambda / cont_ref_nm) ** cont_slope``. The noise is
    Gaussian, with a signal-to-noise ratio ``snr`` per pixel at the
    continuum level (``None``: noiseless). With ``noise_floor = None`` (the
    default) the error is the same at any flux, ``cont / snr``; with a
    value ``a`` between 0 and 1 it follows the flux, as for photon noise
    of the source on a background: sigma = (cont / snr) sqrt(a + (1 - a)
    F / cont), with F the noiseless flux, so that ``a`` is the variance at
    zero flux relative to the continuum.

    ``resol`` is the resolving power, either a single value or ``None`` (no
    instrumental broadening). ``oversample`` is the number of sub-pixels on
    which the model is computed before being averaged over each pixel.

    ``emission_lines`` are added to the continuum (the absorbers absorb
    them too). With ``lyman_limit`` the H I components (those with Lyman
    lines) also absorb the continuum below their Lyman limit (off by
    default, so that the mocks made before it are unchanged). ``n_spikes`` single-pixel positive spikes of
    ``spike_amplitude`` times the local error are added after the noise, at
    random positions; their indices are not part of the truth but are
    reproducible from the seed (see :func:`astrocook.mocks.spike_indices`).
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
    emission_lines: List[EmissionLine] = dataclasses.field(default_factory=list)
    n_spikes: int = 0
    spike_amplitude: float = 10.0
    lyman_limit: bool = False
    noise_floor: Optional[float] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "MockTruth":
        data = dict(data)
        data["components"] = [MockComponent(**c) for c in data["components"]]
        data["emission_lines"] = [EmissionLine(**e)
                                  for e in data.get("emission_lines", [])]
        return cls(**data)

    def to_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_json(cls, path: str) -> "MockTruth":
        with open(path) as f:
            return cls.from_dict(json.load(f))
