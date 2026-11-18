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
from typing import Any, Callable, Dict, List, Tuple

import numpy as np

from astrocook.core.identification import META_KEY as IDENTIFICATION_KEY
from astrocook.core.identification import Identification

logger = logging.getLogger(__name__)

#: The version of the schema of the Bayesian samples
BAYESIAN_SCHEMA_VERSION = 1
BAYESIAN_SCHEMA_NAME = "astrocook.bayesian_samples"


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


#: The structured data of the archive
DATA: Tuple[ArchiveDatum, ...] = (
    ArchiveDatum("identification", "spec", IDENTIFICATION_KEY, "json",
                 _write_identification, _read_identification,
                 "the identification of the regions"),
    ArchiveDatum("bayesian_samples", "systs", "bayesian_results", "npz",
                 _write_bayesian, _read_bayesian, "the samples of the Bayesian fit"),
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
        meta = owner._data.meta if owner is not None else {}
        value = meta.get(datum.meta_key)
        if value is None or (isinstance(value, dict) and not value):
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
    for datum in DATA:
        path = os.path.join(folder, _file_name(base_name, datum, suffix))
        if not os.path.exists(path):
            continue
        try:
            found[datum.owner][datum.meta_key] = datum.read(path)
        except (ValueError, OSError, KeyError, TypeError) as e:
            logger.warning(f"The archive holds {datum.what}, but it cannot be read "
                           f"({e}): the session is opened without it.")
    for owner, extra in ((spec, found["spec"]), (systs, found["systs"])):
        if owner is not None and extra:
            owner._data.meta.update(extra)
