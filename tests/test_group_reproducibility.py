"""The groups of the fits do not depend on the order of a set of UUIDs (F175,
phase 7, PR 7c.4d-ter-bis).

``metals_auto`` gave different numbers in two runs on the same input: in
:meth:`SystemList.get_connected_group` the window of the candidates of the
second line of a multiplet was centred at the redshift of the last member
of the group in the order of a set of random UUIDs (a loop variable reused),
so a component blended with that line was in the group or not by that
order, which changes from run to run (with the UUIDs and the hash seed of
the process).
"""
import json
import os
import subprocess
import sys
import textwrap

from astrocook.core.atomic_data import ATOM_DATA
from astrocook.core.structures import ComponentData, SystemListData
from astrocook.core.system_list import SystemList

C_KMS = 299792.458
Z_SIIV = 2.5


def _components(uuids=None):
    """A SiIV doublet, a Ly-a far from it and a Ly-a blended (5 km/s) with
    SiIV 1402 only (1934 km/s from SiIV 1393, outside the 500 km/s window
    of the candidates of that line)."""
    z_blend = (ATOM_DATA["SiIV_1402"]["wave"] * (1 + Z_SIIV) * (1 + 5.0 / C_KMS)
               / ATOM_DATA["Ly_a"]["wave"] - 1)
    specs = [("SiIV", Z_SIIV), ("Ly_a", 3.2), ("Ly_a", z_blend)]
    comps = []
    for i, (series, z) in enumerate(specs):
        extra = {} if uuids is None else {"uuid": uuids[i]}
        comps.append(ComponentData(id=i + 1, z=z, dz=None, logN=13.0, dlogN=None, b=10.0,
                                   db=None, series=series, **extra))
    return comps


def test_the_second_line_of_a_multiplet_finds_its_blend_whatever_the_order():
    # Before F175 it was found in about half of these lists, by the order of
    # the set of the two seeds
    for _ in range(40):
        systs = SystemList(SystemListData(components=_components()))
        siiv, far, blend = systs.components
        group = systs.get_connected_group([siiv.uuid, far.uuid])
        assert group == {siiv.uuid, far.uuid, blend.uuid}
        assert systs.get_connected_group([siiv.uuid]) == {siiv.uuid, blend.uuid}
        with systs.fitting_context([siiv.uuid, far.uuid], group_depth=0):
            assert systs.constraint_model._active_uuids == {siiv.uuid, far.uuid}
        with systs.fitting_context([siiv.uuid, far.uuid], group_depth=1):
            assert systs.constraint_model._active_uuids == group


def test_fit_groups_break_the_ties_by_the_order_of_the_list():
    # Two components at the same z and of the same series: the seeds and the
    # groups in the order of the list, not of the set
    comps = [ComponentData(id=i + 1, z=Z_SIIV, dz=None, logN=13.0, dlogN=None, b=10.0,
                           db=None, series="SiIV") for i in range(4)]
    for _ in range(10):
        systs = SystemList(SystemListData(components=comps))
        uuids = [c.uuid for c in systs.components]
        assert systs.fit_groups(set(uuids), max_depth=0) == [[u] for u in uuids]
        assert systs.fit_groups(set(uuids), max_depth=1) == [uuids]


_SCRIPT = textwrap.dedent("""
    import json, sys
    sys.path.insert(0, {tests!r})
    from test_group_reproducibility import _components
    from astrocook.core.structures import SystemListData
    from astrocook.core.system_list import SystemList
    uuids = ["u-siiv", "u-far", "u-blend"]
    systs = SystemList(SystemListData(components=_components(uuids)))
    groups = [sorted(systs.get_connected_group(seeds))
              for seeds in (["u-siiv", "u-far"], ["u-far", "u-siiv"])]
    print(json.dumps(groups))
""")


def test_two_processes_with_different_hash_seeds_give_the_same_groups():
    # The order of a set of strings follows the hash seed of the process:
    # with fixed UUIDs, before F175 some of these seeds lost the blend
    script = _SCRIPT.format(tests=os.path.dirname(__file__))
    outputs = set()
    for seed in range(6):
        env = dict(os.environ, PYTHONHASHSEED=str(seed))
        out = subprocess.run([sys.executable, "-c", script], env=env, check=True,
                             capture_output=True, text=True).stdout
        outputs.add(out.strip())
    assert len(outputs) == 1
    (groups,) = [json.loads(o) for o in outputs]
    assert groups == [["u-blend", "u-far", "u-siiv"]] * 2
