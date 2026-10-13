"""Shared pytest configuration and fixtures.

Numerical snapshots
-------------------
The ``snapshot`` fixture compares the output of the current code with
reference arrays stored in ``tests/snapshots/<name>.npz``. These are
characterization tests: they record what the code does now, not what is
correct, and catch unintended numerical changes during refactoring.

Reference files are never written automatically. To create or update them
after an intended change, run::

    pytest --update-snapshots

(add ``--run-slow`` to include the slow snapshots) and review the resulting
diff before committing it.
"""
from pathlib import Path

import numpy as np
import astropy.units as au
import pytest

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"
DATA_DIR = Path(__file__).parent / "data"

# Real spectrum used by the snapshots of continuum, identification and fit.
# The archive does not store the emission redshift: value given by the
# maintainer.
J0103_SESSION = DATA_DIR / "J0103-1305_metallines_lya_new.acs2"
J0103_Z_EM = 4.072
# Lyman-alpha forest, Ly-alpha emission (616.7 nm) and red side, including
# the default rest-frame regions of fit_powerlaw (128-136 nm).
J0103_WINDOW_NM = (570.0, 700.0)
# Resolving power of the J0103-1305 spectrum (given by the maintainer). The
# archive does not store it (saved by an older version of the software).
J0103_RESOL = 45000.0


def pytest_addoption(parser):
    parser.addoption(
        "--update-snapshots", action="store_true", default=False,
        help="Write the numerical snapshots in tests/snapshots/ instead of "
             "comparing against them.")
    parser.addoption(
        "--run-slow", action="store_true", default=False,
        help="Also run the tests marked 'slow' (e.g. the full chain on "
             "J0103-1305, about 15 minutes).")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip = pytest.mark.skip(reason="slow test: run with --run-slow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)


def _to_arrays(values):
    """Convert a dict of arrays/Quantities to plain arrays (units as strings)."""
    arrays = {}
    for key, value in values.items():
        if isinstance(value, au.Quantity):
            arrays[key] = np.asarray(value.value)
            arrays[f"{key}__unit"] = np.array(value.unit.to_string())
        elif value is None:
            arrays[key] = np.array("None")
        else:
            arrays[key] = np.asarray(value)
    return arrays


class Snapshot:
    def __init__(self, update: bool):
        self.update = update

    def check(self, name: str, values: dict, rtol: float = 1e-10, atol: float = 0.0):
        """
        Compare ``values`` with the snapshot ``name``.

        Parameters
        ----------
        name : str
            Snapshot file name (without extension) in ``tests/snapshots/``.
        values : dict
            Arrays, Quantities (value and unit are both checked) or ``None``.
        rtol, atol : float
            Tolerances for floating-point arrays (``np.testing.assert_allclose``).
            NaNs must match exactly.
        """
        path = SNAPSHOT_DIR / f"{name}.npz"
        arrays = _to_arrays(values)

        if self.update:
            SNAPSHOT_DIR.mkdir(exist_ok=True)
            np.savez_compressed(path, **arrays)
            return

        if not path.exists():
            pytest.fail(f"Missing snapshot {path.name}: run pytest --update-snapshots")

        with np.load(path, allow_pickle=False) as ref:
            assert sorted(ref.files) == sorted(arrays), (
                f"Snapshot {name}: keys differ (stored {sorted(ref.files)}, "
                f"got {sorted(arrays)})")
            for key in ref.files:
                expected, actual = ref[key], arrays[key]
                assert expected.shape == actual.shape, (
                    f"Snapshot {name}[{key}]: shape {actual.shape}, "
                    f"expected {expected.shape}")
                if np.issubdtype(expected.dtype, np.inexact):
                    np.testing.assert_allclose(
                        actual, expected, rtol=rtol, atol=atol, equal_nan=True,
                        err_msg=f"Snapshot {name}[{key}]")
                else:
                    np.testing.assert_array_equal(
                        actual, expected, err_msg=f"Snapshot {name}[{key}]")


@pytest.fixture(scope="session")
def j0103_session():
    """The J0103-1305 session (spectrum with continuum, and system list)."""
    from astrocook.core.session import SessionV2
    session = SessionV2.open_new(str(J0103_SESSION), "J0103-1305", None, "auto")
    return session[0] if isinstance(session, (list, tuple)) else session


@pytest.fixture(scope="session")
def j0103_window(j0103_session):
    """
    Raw J0103-1305 spectrum in ``J0103_WINDOW_NM``: x, xmin, xmax, y, dy and
    resol only, with ``z_em = J0103_Z_EM``. Saved continuum and model are
    left out, so that each stage can be run from scratch.
    """
    from astrocook.core.spectrum import SpectrumV2
    from astrocook.core.structures import DataColumnV2, SpectrumDataV2

    spec = j0103_session.spec
    x = spec.x.value
    sel = (x > J0103_WINDOW_NM[0]) & (x < J0103_WINDOW_NM[1])

    def col(q):
        return DataColumnV2(q.value[sel], q.unit)

    data = SpectrumDataV2(
        x=col(spec.x), xmin=col(spec.xmin), xmax=col(spec.xmax),
        y=col(spec.y), dy=col(spec.dy),
        aux_cols={"resol": col(spec.get_column("resol"))},
        z_em=J0103_Z_EM,
    )
    return SpectrumV2(data)


@pytest.fixture(scope="session")
def j0103_window_with_continuum(j0103_session, j0103_window):
    """
    ``j0103_window`` plus the continuum and the absorbed-region mask saved
    in the session (``cont`` and ``abs_mask = ~mask_unabs``). These are the
    frozen inputs of the identification and fit stages: a change in the
    continuum algorithm does not propagate to their snapshots.
    """
    from astrocook.core.spectrum import SpectrumV2
    from astrocook.core.structures import DataColumnV2
    import dataclasses

    spec = j0103_session.spec
    x = spec.x.value
    sel = (x > J0103_WINDOW_NM[0]) & (x < J0103_WINDOW_NM[1])
    cont = spec.get_column("cont")
    unabs = spec.get_column("mask_unabs").value[sel].astype(bool)

    aux = dict(j0103_window._data.aux_cols)
    aux["cont"] = DataColumnV2(cont.value[sel], cont.unit)
    aux["abs_mask"] = DataColumnV2(~unabs, au.dimensionless_unscaled)
    return SpectrumV2(dataclasses.replace(j0103_window._data, aux_cols=aux))


@pytest.fixture(scope="session")
def j0103_window_for_fit(j0103_window_with_continuum):
    """
    ``j0103_window_with_continuum`` with the global resolving power set to
    ``J0103_RESOL`` (as with "Set properties" in the GUI). The all-NaN
    ``resol`` column of the archive is kept, so the fitter uses the global
    value.
    """
    return j0103_window_with_continuum.with_properties(
        z_em=J0103_Z_EM, resol=J0103_RESOL, overwrite_resol_col=False)


@pytest.fixture
def snapshot(request):
    return Snapshot(update=request.config.getoption("--update-snapshots"))
