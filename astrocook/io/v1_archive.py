"""Reading Astrocook V1 session archives (``.acs``) without the V1 code.

An ``.acs`` archive is a tar.gz with one FITS table per structure
(``<name>_spec.fits``, ``<name>_systs.fits``, ...) and the V1 log
(``<name>.json``). The spectrum and the system list are read here directly
into V2 structures. The result is the same as with the former reader,
which went through the V1 classes (``Format.astrocook``, ``Spectrum``,
``SystList``) and converted them to V2, except for what it dropped
(findings register F42, fixed with Guido): the flux unit, and the
``btur``, ``dbtur`` and ``chi2r`` columns of the system table.
The other structures (nodes, lines, models) and the V1 log are not read.

Saving in the ``.acs`` format was removed: Astrocook 1.9.0 can still
write it.
"""
import logging
import os
import tarfile
import tempfile
from typing import Any, Dict, Tuple

import astropy.units as au
import numpy as np
from astropy.io import fits

from astrocook.core.structures import (
    ComponentDataV2, DataColumnV2, SpectrumDataV2, SystemListDataV2,
)

# Columns of the V1 spectrum table read as the core columns
_SPEC_CORE = ['x', 'xmin', 'xmax', 'y', 'dy']
# V1 placed these columns first, in this order, then the others in the
# order of the file ('dx' was skipped)
_SPEC_FIRST = ['cont', 'conv', 'line_mask', 'model', 'deabs', 'resol']
# Flux unit of the V1 spectrum table (as in the V2 'astrocook_fits'
# loader); V1 gave it also to these auxiliary columns
_V1_FLUX_UNIT = au.erg / (au.cm**2 * au.s * au.Angstrom)
_SPEC_FLUX_COLS = ['cont', 'conv', 'model', 'deabs']


class V1ArchiveManager:
    """
    Handles unpacking and cleanup of a V1 Astrocook session (.acs) archive.
    """
    def __init__(self, acs_path: str):
        self.acs_path = acs_path
        self.temp_dir = None

    def unpack(self):
        """Unpacks the .acs (tar.gz) archive into a temporary directory."""
        if not self.acs_path.lower().endswith(('.acs', '.tar.gz')):
            logging.error(f"File is not a recognized V1 archive: {self.acs_path}")
            return None

        try:
            self.temp_dir = tempfile.mkdtemp()
            with tarfile.open(self.acs_path, 'r:gz') as tar:
                # "data" filter: refuse members that would land outside
                # temp_dir (absolute paths, "..", links), F35
                tar.extractall(path=self.temp_dir, filter="data")
            return self.temp_dir
        except Exception as e:
            logging.error(f"Failed to unpack .acs archive {self.acs_path}: {e}")
            self.cleanup()
            return None

    def get_structure_path(self, structure_name: str):
        """Path of the FITS file of a structure (e.g. 'spec'), or None."""
        if not self.temp_dir:
            return None
        # V1 names the files <name>_<structure>.fits
        for fname in os.listdir(self.temp_dir):
            if fname.lower().endswith(f'_{structure_name}.fits'):
                return os.path.join(self.temp_dir, fname)
        return None

    def cleanup(self):
        """Removes the temporary directory created by :meth:`unpack`."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            import shutil
            shutil.rmtree(self.temp_dir)
            logging.debug(f"Cleaned up temporary directory: {self.temp_dir}")


def _v1_meta(hdr: fits.Header) -> fits.Header:
    """Header kept by V1: the first 50 cards, without the long keywords
    that are not Astrocook ones (``AC ...``)."""
    meta = hdr[:50]
    for key in [k for k in meta if 'AC' not in k and len(k) > 8]:
        del meta[key]
    return meta


def read_v1_spectrum(path: str) -> SpectrumDataV2:
    """
    Read the spectrum table of a V1 archive (``<name>_spec.fits``).

    Parameters
    ----------
    path : str
        Path of the FITS file.

    Returns
    -------
    SpectrumDataV2
        x, xmin, xmax in nm; y, dy, cont, conv, model and deabs in
        erg/(cm2 s Angstrom), the other columns without units; ``z_em`` 0
        (V1 did not store it in the table). If the table has no ``resol``
        column, one filled with NaN comes first among the auxiliary
        columns.
    """
    with fits.open(path) as hdul:
        hdr = hdul[1].header
        data = hdul[1].data
        names = list(data.columns.names)

        def core(name, unit, description):
            return DataColumnV2(np.array(data[name], dtype=float), unit,
                                description=description)

        x = core('x', au.nm, "Channels")
        xmin = core('xmin', au.nm, "Lower channel limit")
        xmax = core('xmax', au.nm, "Upper channel limit")
        y = core('y', _V1_FLUX_UNIT, "Flux density")
        dy = core('dy', _V1_FLUX_UNIT, "Error on Flux")

        order = [c for c in _SPEC_FIRST if c in names]
        order += [c for c in names
                  if c not in _SPEC_CORE + ['dx'] and c not in order]
        aux_cols = {}
        if 'resol' not in names:
            aux_cols['resol'] = DataColumnV2(
                np.full(len(x.values), np.nan), au.dimensionless_unscaled,
                description="Spectral Resolution")
        for name in order:
            unit = (_V1_FLUX_UNIT if name in _SPEC_FLUX_COLS
                    else au.dimensionless_unscaled)
            aux_cols[name] = DataColumnV2(
                np.array(data[name]), unit,
                description=f"Auxiliary column: {name}")
        meta = _v1_meta(hdr)

    return SpectrumDataV2(x=x, xmin=xmin, xmax=xmax, y=y, dy=dy,
                          aux_cols=aux_cols, meta=meta, z_em=0.0)


def read_v1_system_list(path: str) -> SystemListDataV2:
    """
    Read the system table of a V1 archive (``<name>_systs.fits``).

    Parameters
    ----------
    path : str
        Path of the FITS file.

    Returns
    -------
    SystemListDataV2
        One component per row, in order of V1 id, with the constraints of
        the ``AC CONSTR`` header keywords. ``chi2`` is the reduced chi2 of
        V1 (``chi2r``); without the ``btur`` and ``dbtur`` columns (early
        V1 archives) ``btur`` is 0 and ``dbtur`` is None, as in V1.
    """
    with fits.open(path) as hdul:
        hdr = dict(hdul[1].header)
        data = hdul[1].data

        # V1 sorted the rows by id; a repeated id kept its first position
        # and the values of its last row
        rows = {}
        for i in np.argsort(data['id'], kind='stable'):
            rows[data['id'][i]] = i

        names = data.columns.names
        components = []
        id_to_uuid = {}
        for v1_id, i in rows.items():
            row = data[i]
            comp = ComponentDataV2(
                id=np.int64(v1_id),
                z=float(row['z']), dz=float(row['dz']),
                logN=float(row['logN']), dlogN=float(row['dlogN']),
                b=float(row['b']), db=float(row['db']),
                btur=float(row['btur']) if 'btur' in names else 0.0,
                dbtur=float(row['dbtur']) if 'dbtur' in names else None,
                func=str(row['func']), series=str(row['series']),
                chi2=float(row['chi2r']) if 'chi2r' in names else None,
                resol=float(row['resol']),
            )
            id_to_uuid[comp.id] = comp.uuid
            components.append(comp)

    return SystemListDataV2(
        components=components,
        meta={},
        v1_header_constraints=hdr,
        parsed_constraints=parse_v1_fits_constraints(hdr),
        v1_id_to_uuid_map=id_to_uuid,
    )


def parse_v1_fits_constraints(hdr: Dict[str, Any]) -> Dict[Tuple[int, str], Dict[str, Any]]:
    """
    Parses V1 constraint keywords (HIERARCH AC CONSTR) from the FITS header
    and returns a clean V2-compatible constraint dictionary.

    Returns: {(comp_id, param_name): {'is_free': bool, 'expression': str}}
    """
    # Keys will store the final parsed constraints
    parsed_constraints = {}

    # Storage for intermediate parsing (key is the progressive integer [n])
    constraint_data = {}
    # Iterate over the header items to find matching keys
    for key, value in hdr.items():
        if key.startswith('AC CONSTR'):
            try:
                parts = key.split()
                # Assuming key structure is AC CONSTR ID [n]
                # The index [n] is the 4th element (index 3) of the split list
                # The type is the 3rd element (index 2), e.g., 'ID', 'PAR', 'VAL'

                # Check for ID, PAR, and VAL keywords using the progressive index [n]
                if len(parts) >= 4:
                    index_type = parts[2] # 'ID', 'PAR', or 'VAL'
                    n = int(parts[3])      # The progressive integer index [n]

                    if n not in constraint_data:
                        constraint_data[n] = {}

                    if index_type == 'ID':
                        constraint_data[n]['id'] = int(value)
                    elif index_type == 'PAR':
                        constraint_data[n]['par'] = value.lower() # Normalize parameter name
                    elif index_type == 'VAL':
                        constraint_data[n]['val'] = value

            except (ValueError, IndexError):
                # Ignore malformed or unexpected header keywords
                continue

    # --- Convert Intermediate Data to Final V2 Structure ---
    for n, data in constraint_data.items():
        if 'id' in data and 'par' in data and 'val' in data:
            comp_id = data['id']
            param_name = data['par']
            val = data['val']

            # 1. Determine Freezing (vary=False)
            is_free = True
            if val == '' or val is None:
                is_free = False # VAL is empty string means frozen

            # 2. Determine Expression (link/fixed value)
            expression = val if val != '' else None

            # V2 Key: (comp_id, param_name)
            parsed_constraints[(comp_id, param_name)] = {
                'is_free': is_free,
                'expression': expression
            }

    logging.info(f"Parsed {len(parsed_constraints)} constraints from FITS header metadata.")
    return parsed_constraints
