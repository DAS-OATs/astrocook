"""Full analysis chain, run in a separate process by the full-chain snapshots.

The chain is the one a user runs from the GUI menus, with the recipe
defaults: Edit > Set properties (z_em, R), Continuum > Auto-estimate,
Absorbers > Auto-detect doublets (red side), Absorbers > Auto-fit Ly-alpha
forest.

The result of the automatic pipelines depends on the order in which
``refit_all`` picks its seeds from a set of component UUIDs, i.e. on the
random UUIDs and on the string hash seed (findings register F19,
``ROADMAP.md``). To get a reproducible result the chain runs in a child
process with ``PYTHONHASHSEED=0`` and with ``uuid.uuid4`` replaced by a
seeded generator. Usage::

    PYTHONHASHSEED=0 python tests/full_chain.py <case> <output.npz>

with ``<case>`` one of :data:`CASES`.
"""
import dataclasses
import logging
import random
import sys
import uuid
from pathlib import Path

import numpy as np

TESTS_DIR = Path(__file__).parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from conftest import J0103_RESOL, J0103_SESSION, J0103_Z_EM  # noqa: E402

# Wider than the window of the stage snapshots: down to 500 nm, to include
# the Lyman-beta emission (520.2 nm) and part of the Lyman-beta forest.
J0103_CHAIN_WINDOW_NM = (500.0, 700.0)

# Mock: z_em = 2.5, Lyman-beta emission at 359 nm, Lyman-alpha at 425.5 nm.
MOCK_Z_EM = 2.5
MOCK_WINDOW_NM = (355.0, 460.0)

#: Continuum stored every STEP pixels (it is smooth)
STEP = 10

CASES = ("j0103", "mock")


def seed_uuids(seed: int = 0) -> None:
    """Replace ``uuid.uuid4`` with a reproducible generator."""
    rng = random.Random(seed)
    uuid.uuid4 = lambda: uuid.UUID(int=rng.getrandbits(128), version=4)


def j0103_session():
    """Raw J0103-1305 spectrum in the chain window; z_em and R are set by
    the chain itself, as a user would do."""
    from astrocook.core.session import SessionV2
    from astrocook.core.spectrum import SpectrumV2
    from astrocook.core.structures import DataColumnV2, SpectrumDataV2

    loaded = SessionV2.open_new(str(J0103_SESSION), "J0103-1305", None, "auto")
    spec = loaded.spec
    x = spec.x.value
    sel = (x > J0103_CHAIN_WINDOW_NM[0]) & (x < J0103_CHAIN_WINDOW_NM[1])

    def col(q):
        return DataColumnV2(q.value[sel], q.unit)

    data = SpectrumDataV2(
        x=col(spec.x), xmin=col(spec.xmin), xmax=col(spec.xmax),
        y=col(spec.y), dy=col(spec.dy),
        aux_cols={"resol": col(spec.get_column("resol"))})
    session = SessionV2("J0103-1305", None, spec=SpectrumV2(data))
    return session.edit.set_properties(z_em=str(J0103_Z_EM),
                                       resol=str(J0103_RESOL))


def mock_truth():
    """Lyman-alpha + Lyman-beta forest, three CIV systems on the red side
    (one with two components) and a SiIV doublet in the forest, on a
    power-law continuum with emission lines."""
    from astrocook.mocks import MockComponent, MockTruth
    from astrocook.mocks.populations import lya_forest
    from astrocook.mocks.truth import QSO_EMISSION_LINES

    rng = np.random.default_rng(7)
    z_min = MOCK_WINDOW_NM[0] / 121.567 - 1.0
    forest = [dataclasses.replace(c, series="Ly_ab")
              for c in lya_forest(rng, z_min, MOCK_Z_EM)]
    metals = [MockComponent("CIV", 1.7600, 13.8, 10.0),
              MockComponent("CIV", 1.8300, 13.4, 8.0),
              MockComponent("CIV", 1.83012, 13.9, 12.0),
              MockComponent("CIV", 1.9000, 14.2, 15.0),
              MockComponent("SiIV", 1.9000, 13.3, 12.0)]
    return MockTruth(components=forest + metals,
                     lam_min_nm=MOCK_WINDOW_NM[0], lam_max_nm=MOCK_WINDOW_NM[1],
                     dv_kms=2.5, resol=45000.0, snr=30.0, z_em=MOCK_Z_EM,
                     seed=7, emission_lines=QSO_EMISSION_LINES, cont_slope=-1.0)


def mock_session():
    """Mock spectrum; z_em and R are already set by the generator."""
    from astrocook.core.session import SessionV2
    from astrocook.mocks import generate
    return SessionV2("mock", None, spec=generate(mock_truth()))


def run_chain(session):
    """Continuum, doublets on the red side, Lyman-alpha forest."""
    continuum = session.continuum.estimate_auto()
    doublets = continuum.absorbers.doublets_auto()
    forest = doublets.absorbers.lya_auto()
    return continuum, doublets, forest


def _components(systs, prefix):
    comps = sorted(systs.components, key=lambda c: (c.z, c.series))
    out = {f"{prefix}_{key}": np.array([getattr(c, key) for c in comps], dtype=float)
           for key in ("z", "dz", "logN", "dlogN", "b", "db")}
    out[f"{prefix}_series"] = np.array([c.series for c in comps])
    return out


def chain_arrays(continuum, doublets, forest):
    """Arrays recorded by the snapshots, one group per stage."""
    spec = continuum.spec
    cont = spec.get_column("cont").value
    out = {"cont": cont[::STEP],
           "abs_mask": spec.get_column("abs_mask").value.astype(bool)}
    out.update(_components(doublets.systs, "doublets"))
    out.update(_components(forest.systs, "forest"))
    model = forest.spec.get_column("model").value
    out["model_norm"] = (model / forest.spec.get_column("cont").value).astype(np.float32)
    return out


def main(case: str, path: str) -> None:
    logging.basicConfig(level=logging.WARNING)
    seed_uuids(0)
    session = {"j0103": j0103_session, "mock": mock_session}[case]()
    np.savez(path, **chain_arrays(*run_chain(session)))


if __name__ == "__main__":
    main(*sys.argv[1:])
