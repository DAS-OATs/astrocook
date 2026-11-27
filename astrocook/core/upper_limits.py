"""The upper limits on the column of a species not detected (phase 7, PR
7d.7; findings register F156).

A limit says that a species at a redshift, if there, has a column below
``logN_lim`` at a confidence of ``sigma``: it is not a component, and the
model never includes it. The limits live in the metadata of the system list
under :data:`META_KEY` (decided with Guido, 27 November 2026: the System
Inspector shows the system list, and the identification is rewritten by each
run of the identification), so that they follow the history of the session
and the undo with no extra code; in the archive they are a file of their
own, ``<name>_upper_limits.json`` (``astrocook.io.archive_data``).

A limit describes the data and the model on its pixels when it was computed:
it is removed when they change (:mod:`astrocook.fitting.upper_limits`,
:func:`~astrocook.fitting.upper_limits.drop_stale`, Guido, 27 November
2026), not recomputed. Its ``components`` and ``data`` are the fingerprints
of that state.

Computed by :func:`astrocook.fitting.upper_limits.upper_limit`, through the
recipe ``absorbers.upper_limit``.
"""
from __future__ import annotations

import dataclasses
import uuid as _uuid
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

#: The key of the limits in the metadata of the system list
META_KEY = "upper_limits"
#: The version of the schema of the file of the archive
SCHEMA_VERSION = 1
#: The name of the schema, written in the file
SCHEMA_NAME = "astrocook.upper_limits"

#: The methods of a limit
METHODS = ("profile", "ew")


@dataclass(frozen=True)
class UpperLimit:
    """An upper limit on the column of a species at a redshift."""
    series: str
    """ The species (all its lines) or a transition, as a component's. """
    z: float
    """ The redshift. """
    b: float
    """ The Doppler parameter assumed, fixed (km/s). """
    logN_lim: float
    """ The limit, log cm^-2. """
    sigma: float
    """ The confidence in sigma: the logN where the chi2 rises by sigma^2
    over its minimum (``profile``), or sigma times the error of the column
    from the equivalent width (``ew``). """
    method: str
    """ ``profile`` or ``ew``. """
    transitions: Tuple[str, ...]
    """ The lines that constrain it (in the spectrum, with valid pixels). """
    skipped: Tuple[str, ...] = ()
    """ The lines of ``series`` left out: outside the spectrum, or without
    valid pixels. """
    ranges_nm: Tuple[Tuple[float, float], ...] = ()
    """ The pixels of the limit, as wavelength ranges (nm). """
    n_pix: int = 0
    """ The pixels of the limit. """
    b_inst: float = 0.0
    """ The instrumental b at the strongest line (km/s; 0 if unknown). """
    logN_best: Optional[float] = None
    """ The logN at the minimum of the chi2 (``profile``), or the column of
    the equivalent width measured (``ew``); None below the scan. """
    significance: float = 0.0
    """ How much better the best logN is than no absorption, in sigma
    (sqrt of the difference of the chi2; for ``ew``, the measured column
    over its error). Above ``sigma`` the species is detected and the
    limit is not meaningful (the recipe warns). """
    refitted: Tuple[str, ...] = ()
    """ The uuids of the components refitted at each step (``profile``):
    those with a line on the pixels. """
    ref_uuid: Optional[str] = None
    """ The component whose redshift the limit took, if any. """
    reach_nm: Tuple[Tuple[float, float], ...] = ()
    """ The ranges (nm) where a component changes the model on the pixels:
    those of the pixels widened by the instrumental kernel; the components
    with a line there are those of the fingerprint. """
    components: str = ""
    """ The fingerprint of the components with a line in ``reach_nm``
    (their parameters and constraints) when the limit was computed. """
    data: str = ""
    """ The fingerprint of the data on the pixels (wavelength, flux, error,
    continuum and resolution) when the limit was computed. """
    recipe: Optional[str] = None
    """ The recipe that computed it. """
    time: str = ""
    """ When, ISO 8601 in UTC. """
    uuid: str = field(default_factory=lambda: str(_uuid.uuid4()))
    """ Its identifier (``absorbers.delete_upper_limit``). """

    def to_dict(self) -> Dict[str, Any]:
        out = {f: getattr(self, f) for f in self.__dataclass_fields__}
        for key in ("transitions", "skipped", "refitted"):
            out[key] = list(out[key])
        for key in ("ranges_nm", "reach_nm"):
            out[key] = [list(r) for r in out[key]]
        return out

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> 'UpperLimit':
        """The limit of :meth:`to_dict`; a key that this version does not
        know is ignored."""
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        for key in ("transitions", "skipped", "refitted"):
            known[key] = tuple(str(v) for v in known.get(key, ()))
        for key in ("ranges_nm", "reach_nm"):
            known[key] = tuple((float(lo), float(hi)) for lo, hi in known.get(key, ()))
        return cls(**known)

    def label(self) -> str:
        """``logN < 12.34`` (the confidence apart)."""
        return f"< {self.logN_lim:.2f}"


def limits_of(systs: Any) -> Tuple[UpperLimit, ...]:
    """The limits of a system list (none if it has none, or no list)."""
    if systs is None:
        return ()
    return tuple(systs._data.meta.get(META_KEY, ()))


def with_limits(systs: Any, limits: Tuple[UpperLimit, ...]) -> Any:
    """A new system list with ``limits`` (the components unchanged)."""
    from astrocook.core.system_list import SystemList
    meta = dict(systs._data.meta)
    meta[META_KEY] = tuple(limits)
    return SystemList(dataclasses.replace(systs._data, meta=meta))


def to_file_dict(limits: Tuple[UpperLimit, ...]) -> Dict[str, Any]:
    """The content of the file of the archive."""
    return {"schema": SCHEMA_NAME, "schema_version": SCHEMA_VERSION,
            "limits": [lim.to_dict() for lim in limits]}


def from_file_dict(d: Dict[str, Any]) -> Tuple[UpperLimit, ...]:
    """The limits of :func:`to_file_dict`.

    Raises
    ------
    ValueError
        If it is not a file of limits, or a later version of its schema.
    """
    if d.get("schema") != SCHEMA_NAME:
        raise ValueError(f"not the upper limits: schema {d.get('schema')!r}")
    version = int(d.get("schema_version", 0))
    if version > SCHEMA_VERSION:
        raise ValueError(f"the version {version} of the schema is not known to this "
                         f"Astrocook (it reads up to {SCHEMA_VERSION})")
    return tuple(UpperLimit.from_dict(x) for x in d.get("limits", ()))
