"""The structured data of an ``.acs2`` archive (phase 7, PR 7c.3; F28).

Data that a FITS header cannot hold live in the metadata of their owner in
memory (the spectrum or the system list), so that they follow the history
of the session and the undo with no extra code, and in the archive as
files of their own in the tar, one per datum, named
``<name>_<datum>.<ext>`` (``<name>_<datum>_orig.<ext>`` for the original
data of the session), each with the name and the version of its schema:

- ``identification`` (``.json``): the identification of the absorbed
  regions, :class:`~astrocook.core.identification.Identification`, in the
  metadata of the spectrum;
- ``bayesian_samples`` (``.npz``, without pickle): the samples of
  ``fit_bayesian`` (``bayesian_results`` in the metadata of the system
  list, read by the corner plot).
- ``fit_records`` (``.json``): the statistics of the group fits (phase 7,
  PR 7d.9), :class:`~astrocook.core.structures.FitRecord`, each written
  once, with the record of each component; in memory they are not in the
  metadata but in the components (``ComponentData.fit``), which
  :attr:`ArchiveDatum.get` and :attr:`ArchiveDatum.put` read and fill.
- ``upper_limits`` (``.json``): the upper limits on the column of the
  species not detected (phase 7, PR 7d.7),
  :class:`~astrocook.core.upper_limits.UpperLimit`, in the metadata of the
  system list.

A datum is added by an entry of :data:`DATA`. The archives saved before
this module open unchanged (they have no such files), and a version that
does not know a file ignores it: the loader of 2.7 looks only for the
``_spec.fits``, ``_systs.fits`` and ``_meta.json`` members (checked with
2.7.1 on an archive of this version). A file that cannot be read (a later
version of its schema, a damaged file) is left out with a WARNING: the
session opens without it.
"""
import json
import logging
import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from astrocook.core.identification import META_KEY as IDENTIFICATION_KEY
from astrocook.core.identification import Identification
from astrocook.core.structures import FitRecord
from astrocook.core.upper_limits import META_KEY as UPPER_LIMITS_KEY
from astrocook.core.upper_limits import from_file_dict as _limits_from_dict
from astrocook.core.upper_limits import to_file_dict as _limits_to_dict

logger = logging.getLogger(__name__)

#: The version of the schema of the Bayesian samples
BAYESIAN_SCHEMA_VERSION = 1
BAYESIAN_SCHEMA_NAME = "astrocook.bayesian_samples"
#: The version of the schema of the statistics of the fits
FIT_RECORDS_SCHEMA_VERSION = 1
FIT_RECORDS_SCHEMA_NAME = "astrocook.fit_records"


@dataclass(frozen=True)
class ArchiveDatum:
    """A structured datum of the archive."""
    name: str
    """Its part of the file name, ``<name>_<datum>.<ext>``."""
    owner: str
    """``"spec"`` or ``"systs"``: whose metadata hold it."""
    meta_key: str
    ext: str
    write: Callable[[Any, str], None]
    """Writes the value to a path."""
    read: Callable[[str], Any]
    """Reads the value from a path; raises ``ValueError`` (or ``OSError``,
    ``KeyError``) if it cannot."""
    what: str
    """What it is, in words, for the messages."""
    get: Optional[Callable[[Any], Any]] = None
    """Reads the value from its owner; None: ``meta[meta_key]``."""
    put: Optional[Callable[[Any, Any], None]] = None
    """Gives a value read to its owner, just read from the archive (in
    place: no one else holds it yet); None: into ``meta[meta_key]``."""


# --- identification (JSON) ---

def _write_identification(value: Any, path: str) -> None:
    if not isinstance(value, Identification):
        raise TypeError(f"not an Identification: {type(value).__name__}")
    with open(path, "w") as f:
        json.dump(value.to_dict(), f, indent=1)


def _read_identification(path: str) -> Identification:
    with open(path) as f:
        return Identification.from_dict(json.load(f))


# --- Bayesian samples (npz) ---

def _write_bayesian(value: Any, path: str) -> None:
    """``{uuid: {'samples', 'labels', 'engine'}}``; the components sampled
    together share one entry (``fit_bayesian``), written once."""
    groups: Dict[int, Tuple[dict, List[str]]] = {}
    skipped = 0
    for uuid, res in dict(value).items():
        if not isinstance(res, dict) or "samples" not in res:
            skipped += 1
            continue
        groups.setdefault(id(res), (res, []))[1].append(str(uuid))
    if skipped:
        logger.warning(f"{skipped} entries of the Bayesian samples are not in the form "
                       "written by fit_bayesian: they are not saved.")
    arrays: Dict[str, np.ndarray] = {
        "schema": np.array(BAYESIAN_SCHEMA_NAME),
        "schema_version": np.array(BAYESIAN_SCHEMA_VERSION)}
    for i, (res, uuids) in enumerate(groups.values()):
        samples = np.asarray(res["samples"], dtype=float)
        labels = list(res.get("labels") or [f"p_{j}" for j in range(samples.shape[-1])])
        arrays[f"samples_{i}"] = samples
        arrays[f"labels_{i}"] = np.array(labels, dtype=str)
        arrays[f"engine_{i}"] = np.array(str(res.get("engine", "")))
        arrays[f"uuids_{i}"] = np.array(uuids, dtype=str)
    with open(path, "wb") as f:
        np.savez_compressed(f, **arrays)


def _read_bayesian(path: str) -> Dict[str, dict]:
    with np.load(path, allow_pickle=False) as npz:
        if str(npz["schema"]) != BAYESIAN_SCHEMA_NAME:
            raise ValueError(f"not Bayesian samples: schema {str(npz['schema'])!r}")
        version = int(npz["schema_version"])
        if version > BAYESIAN_SCHEMA_VERSION:
            raise ValueError(f"the version {version} of the schema is not known to this "
                             f"Astrocook (it reads up to {BAYESIAN_SCHEMA_VERSION})")
        out: Dict[str, dict] = {}
        i = 0
        while f"samples_{i}" in npz.files:
            res = {"samples": npz[f"samples_{i}"],
                   "labels": [str(s) for s in npz[f"labels_{i}"]],
                   "engine": str(npz[f"engine_{i}"])}
            for uuid in npz[f"uuids_{i}"]:
                out[str(uuid)] = res
            i += 1
    return out


# --- statistics of the group fits (JSON; 7d.9) ---

def _get_fit_records(systs: Any) -> Any:
    """``{'records': {id: FitRecord}, 'components': {uuid: id}}``; None
    if no component holds a record."""
    components = {c.uuid: c.fit.id for c in systs.components if c.fit is not None}
    if not components:
        return None
    records = {}
    for c in systs.components:
        if c.fit is not None:
            records.setdefault(c.fit.id, c.fit)
    return {"records": records, "components": components}


def _write_fit_records(value: Any, path: str) -> None:
    with open(path, "w") as f:
        json.dump({"schema": FIT_RECORDS_SCHEMA_NAME,
                   "schema_version": FIT_RECORDS_SCHEMA_VERSION,
                   "records": [r.to_dict() for r in value["records"].values()],
                   "components": value["components"]}, f, indent=1)


def _read_fit_records(path: str) -> Dict[str, Any]:
    with open(path) as f:
        d = json.load(f)
    if d.get("schema") != FIT_RECORDS_SCHEMA_NAME:
        raise ValueError(f"not the statistics of the fits: schema {d.get('schema')!r}")
    version = int(d.get("schema_version", 0))
    if version > FIT_RECORDS_SCHEMA_VERSION:
        raise ValueError(f"the version {version} of the schema is not known to this "
                         f"Astrocook (it reads up to {FIT_RECORDS_SCHEMA_VERSION})")
    records = {r["id"]: FitRecord.from_dict(r) for r in d["records"]}
    return {"records": records, "components": dict(d["components"])}


def _put_fit_records(systs: Any, value: Dict[str, Any]) -> None:
    """Give each component of ``systs``, just read, its record (frozen
    dataclasses, set in place as the loader sets their uuid)."""
    records, of = value["records"], value["components"]
    for c in systs.components:
        record = records.get(of.get(c.uuid))
        if record is not None:
            object.__setattr__(c, "fit", record)


# --- upper limits (JSON) ---

def _write_upper_limits(value: Any, path: str) -> None:
    with open(path, "w") as f:
        json.dump(_limits_to_dict(tuple(value)), f, indent=1)


def _read_upper_limits(path: str) -> Any:
    with open(path) as f:
        return _limits_from_dict(json.load(f))


#: The structured data of the archive
DATA: Tuple[ArchiveDatum, ...] = (
    ArchiveDatum("identification", "spec", IDENTIFICATION_KEY, "json",
                 _write_identification, _read_identification,
                 "the identification of the regions"),
    ArchiveDatum("bayesian_samples", "systs", "bayesian_results", "npz",
                 _write_bayesian, _read_bayesian, "the samples of the Bayesian fit"),
    ArchiveDatum("fit_records", "systs", "fit_records", "json",
                 _write_fit_records, _read_fit_records, "the statistics of the fits",
                 get=_get_fit_records, put=_put_fit_records),
    ArchiveDatum("upper_limits", "systs", UPPER_LIMITS_KEY, "json",
                 _write_upper_limits, _read_upper_limits, "the upper limits"),
)


def meta_keys(owner: str) -> set:
    """The metadata keys of ``owner`` saved as files: not in the FITS header."""
    return {d.meta_key for d in DATA if d.owner == owner}


def _file_name(base_name: str, datum: ArchiveDatum, suffix: str) -> str:
    return f"{base_name}_{datum.name}{suffix}.{datum.ext}"


def write_data(session: Any, folder: str, base_name: str,
               suffix: str = "") -> List[str]:
    """Write the structured data of ``session`` to ``folder``.

    Returns
    -------
    list of str
        The names of the files written (members of the archive).
    """
    names = []
    owners = {"spec": session.spec, "systs": session.systs}
    for datum in DATA:
        owner = owners[datum.owner]
        if datum.get is not None:
            value = datum.get(owner) if owner is not None else None
        else:
            meta = owner._data.meta if owner is not None else {}
            value = meta.get(datum.meta_key)
        if value is None or (isinstance(value, (dict, tuple, list)) and not value):
            continue
        name = _file_name(base_name, datum, suffix)
        datum.write(value, os.path.join(folder, name))
        names.append(name)
    return names


def read_data(archive_root: str, spec: Any, systs: Any, suffix: str = "") -> None:
    """Add to the metadata of the spectrum and the system list just read
    from an archive (in place: no one else holds them yet) the structured
    data found next to ``archive_root`` (the path of the members without
    ``_spec.fits``).

    A file that cannot be read is left out with a WARNING.
    """
    folder, base_name = os.path.split(archive_root)
    found: Dict[str, Dict[str, Any]] = {"spec": {}, "systs": {}}
    owners = {"spec": spec, "systs": systs}
    for datum in DATA:
        path = os.path.join(folder, _file_name(base_name, datum, suffix))
        if not os.path.exists(path):
            continue
        try:
            value = datum.read(path)
            if datum.put is not None:
                if owners[datum.owner] is not None:
                    datum.put(owners[datum.owner], value)
            else:
                found[datum.owner][datum.meta_key] = value
        except (ValueError, OSError, KeyError, TypeError) as e:
            logger.warning(f"The archive holds {datum.what}, but it cannot be read "
                           f"({e}): the session is opened without it.")
    for owner, extra in ((spec, found["spec"]), (systs, found["systs"])):
        if owner is not None and extra:
            owner._data.meta.update(extra)
