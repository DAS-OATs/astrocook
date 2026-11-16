"""Ground truth of a mock spectrum."""
import dataclasses
import json
from typing import List, Optional, Tuple


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
    shift_kms : float, optional
        Velocity shift of the peak from the systemic redshift ``z_em``
        (km/s, negative to the blue, e.g. the CIV blueshift).
    asymmetry : float, optional
        Between -1 and 1: the half of the line on the blue side of the peak
        has a FWHM of ``fwhm_kms * (1 + asymmetry)``, the red one of
        ``fwhm_kms * (1 - asymmetry)`` (positive: a broader blue wing).
    """
    rest_nm: float
    amplitude: float
    fwhm_kms: float
    shift_kms: float = 0.0
    asymmetry: float = 0.0


#: Lyman-alpha + NV and CIV emission, indicative values
QSO_EMISSION_LINES = [EmissionLine(121.567, 2.0, 4000.0),
                      EmissionLine(123.9, 0.4, 3000.0),
                      EmissionLine(154.9, 0.8, 5000.0)]

#: A richer set of quasar emission lines (phase 7, PR 7a.1), with the
#: shift of CIV to the blue and its broader blue wing, at the systemic
#: z_em: OVI and Ly-b, Ly-a, NV, SiIV+OIV], CIV, CIII]. Indicative values,
#: after the composite spectrum of Vanden Berk et al. (2001, AJ 122, 549)
#: for the positions and the relative strengths, and Richards et al.
#: (2011, AJ 141, 167) for the CIV blueshift (a few hundred to a few
#: thousand km/s). Not calibrated: the mocks of the robustness grid draw
#: their amplitudes, shifts and asymmetries around these values.
QSO_EMISSION_LINES_FULL = [
    EmissionLine(102.572, 0.15, 4000.0),                  # Ly-b
    EmissionLine(103.383, 0.3, 5000.0),                   # OVI 1032, 1038
    EmissionLine(121.567, 2.0, 4000.0),                   # Ly-a
    EmissionLine(124.08, 0.4, 3000.0),                    # NV 1239, 1243
    EmissionLine(139.8, 0.25, 5000.0, shift_kms=-300.0),  # SiIV + OIV]
    EmissionLine(154.906, 0.8, 5000.0, shift_kms=-800.0,
                 asymmetry=0.2),                          # CIV 1548, 1551
    EmissionLine(190.873, 0.4, 5000.0, shift_kms=-200.0)]  # CIII]


@dataclasses.dataclass(frozen=True)
class BALTrough:
    """A broad absorption line (BAL) trough (phase 7, PR 7a.1).

    The optical depth is smooth in the outflow velocity v (km/s, positive
    to the blue of the systemic z_em, relativistic Doppler formula): flat
    at ``tau_max`` between ``v_min_kms + edge_kms`` and ``v_max_kms -
    edge_kms``, falling to zero at ``v_min_kms`` and ``v_max_kms`` as a
    half cosine. Each line of each ion has the optical depth ``tau_max *
    ratio * (f lambda) / (f lambda)_CIV_1548``, as for the same column
    density of every ion (``column_ratios`` all 1, the default) or with
    the given ratios of the columns to that of CIV.

    Parameters
    ----------
    v_min_kms, v_max_kms : float
        Velocity range of the trough (km/s, outflow positive).
    tau_max : float
        Optical depth of CIV 1548 on the flat part.
    edge_kms : float, optional
        Width of each of the two edges (at most half the range).
    ions : tuple of str, optional
        Multiplets of ``STANDARD_MULTIPLETS`` or transitions.
    column_ratios : tuple of float, optional
        Column of each ion relative to CIV; ``None``: all 1.
    """
    v_min_kms: float
    v_max_kms: float
    tau_max: float
    edge_kms: float = 1000.0
    ions: Tuple[str, ...] = ("CIV", "SiIV", "NV", "Ly_a")
    column_ratios: Optional[Tuple[float, ...]] = None


@dataclasses.dataclass(frozen=True)
class TelluricModel:
    """Telluric absorption of a mock (phase 7, PR 7a.1).

    The molecular transmission of the template
    (:mod:`astrocook.mocks.telluric`, SkyCalc, Paranal, airmass 1.1, PWV
    2.5 mm), with its optical depth split by wavelength into an O2 part
    (the gamma, B and A bands) and an H2O part (the rest): the first
    scaled by ``airmass / 1.1``, the second by ``airmass / 1.1 * pwv_mm /
    2.5``, both shifted by ``shift_kms``. Applied in the observed frame,
    before the instrumental convolution.
    """
    airmass: float = 1.1
    pwv_mm: float = 2.5
    shift_kms: float = 0.0


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

    Added in phase 7 (PR 7a.1), all off by default, so that the mocks made
    before them are unchanged: ``bal_troughs`` (:class:`BALTrough`, on the
    systemic ``z_em``), ``telluric`` (:class:`TelluricModel`) and
    ``proximity_kms``, the velocity from ``z_em`` at which the ionizing
    flux of the quasar equals the background in the forest of
    :func:`~astrocook.mocks.lya_forest`: a record of how the forest was
    drawn, for the metrics, not used by :func:`~astrocook.mocks.generate`.
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
    bal_troughs: List[BALTrough] = dataclasses.field(default_factory=list)
    telluric: Optional[TelluricModel] = None
    proximity_kms: Optional[float] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "MockTruth":
        data = dict(data)
        data["components"] = [MockComponent(**c) for c in data["components"]]
        data["emission_lines"] = [EmissionLine(**e)
                                  for e in data.get("emission_lines", [])]
        data["bal_troughs"] = [_bal_from_dict(t) for t in data.get("bal_troughs", [])]
        if data.get("telluric") is not None:
            data["telluric"] = TelluricModel(**data["telluric"])
        return cls(**data)

    def to_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_json(cls, path: str) -> "MockTruth":
        with open(path) as f:
            return cls.from_dict(json.load(f))


def _bal_from_dict(data: dict) -> BALTrough:
    data = dict(data)
    data["ions"] = tuple(data.get("ions", BALTrough.ions))
    if data.get("column_ratios") is not None:
        data["column_ratios"] = tuple(data["column_ratios"])
    return BALTrough(**data)
