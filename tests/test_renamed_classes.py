"""The classes that lost the ``V2`` suffix (phase 6-quater, PR 6): the old
names still import, give the same class, and warn (``astrocook/_renamed.py``)."""
import importlib
import warnings

import pytest

from astrocook._renamed import renamed

# (module, old name, new name): the modules that defined or exported a class
RENAMED = [
    ("astrocook.core.structures", "SessionMetadataV2", "SessionMetadata"),
    ("astrocook.core.structures", "DataColumnV2", "DataColumn"),
    ("astrocook.core.structures", "SpectrumDataV2", "SpectrumData"),
    ("astrocook.core.structures", "ComponentDataV2", "ComponentData"),
    ("astrocook.core.structures", "ParameterConstraintV2", "ParameterConstraint"),
    ("astrocook.core.structures", "SystemListDataV2", "SystemListData"),
    ("astrocook.core.structures", "LogEntryV2", "LogEntry"),
    ("astrocook.core.structures", "LogNoteV2", "LogNote"),
    ("astrocook.core.structures", "HistoryLogV2", "HistoryLog"),
    ("astrocook.core.spectrum", "TableAdapterV2", "TableAdapter"),
    ("astrocook.core.spectrum", "SpectrumV2", "Spectrum"),
    ("astrocook.core.system_list", "SystemListV2", "SystemList"),
    ("astrocook.fitting.voigt_model", "VoigtModelConstraintV2", "VoigtModelConstraint"),
    ("astrocook.fitting.voigt_fitter", "VoigtFitterV2", "VoigtFitter"),
    ("astrocook.recipes.file", "RecipeFileV2", "RecipeFile"),
    ("astrocook.recipes.flux", "RecipeFluxV2", "RecipeFlux"),
    ("astrocook.recipes.continuum", "RecipeContinuumV2", "RecipeContinuum"),
    ("astrocook.recipes.edit", "RecipeEditV2", "RecipeEdit"),
    ("astrocook.recipes.absorbers", "RecipeAbsorbersV2", "RecipeAbsorbers"),
    ("astrocook.gui.main_window", "MainWindowV2", "MainWindow"),
    ("astrocook.session.session", "SessionV2", "Session"),
    # where the old names were also imported from
    ("astrocook", "SessionV2", "Session"),
    ("astrocook.session", "SessionV2", "Session"),
    ("astrocook.core", "SessionV2", "Session"),
    ("astrocook.core", "SpectrumDataV2", "SpectrumData"),
    ("astrocook.core", "SystemListDataV2", "SystemListData"),
    ("astrocook.core.session", "SessionV2", "Session"),  # path of PR 5
]


@pytest.mark.parametrize("module,old,new", RENAMED)
def test_old_name_is_the_new_class_and_warns(module, old, new):
    mod = importlib.import_module(module)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cls = getattr(mod, new)
    assert cls.__name__ == new
    with pytest.warns(DeprecationWarning, match=f"{old}.*{new}"):
        assert getattr(mod, old) is cls


def test_old_name_by_from_import_warns_at_the_importer():
    with pytest.warns(DeprecationWarning) as record:
        from astrocook.core.spectrum import SpectrumV2  # noqa: F401
    assert record[0].filename == __file__


def test_isinstance_with_the_old_name():
    from astrocook.core.spectrum import Spectrum
    from astrocook.core.structures import DataColumn, SpectrumData
    from astropy import units as au
    import numpy as np
    col = DataColumn(np.arange(3.0), au.nm)
    flux = DataColumn(np.ones(3), au.dimensionless_unscaled)
    spec = Spectrum(data=SpectrumData(x=col, xmin=col, xmax=col, y=flux, dy=flux))
    with pytest.warns(DeprecationWarning):
        from astrocook.core.spectrum import SpectrumV2
    assert isinstance(spec, SpectrumV2)


def test_unknown_name_is_an_attribute_error():
    import astrocook.core.spectrum as mod
    with pytest.raises(AttributeError):
        mod.NoSuchThing  # noqa: B018
    with pytest.raises(AttributeError):
        mod.Spectrum2  # noqa: B018


def test_helper_looks_the_new_name_up_when_asked():
    import sys
    import types
    mod = types.ModuleType("fake_renamed_module")
    sys.modules[mod.__name__] = mod
    try:
        mod.__getattr__ = renamed(mod.__name__, {"OldV2": "New"})  # type: ignore[attr-defined]
        mod.New = object()  # type: ignore[attr-defined]
        with pytest.warns(DeprecationWarning):
            assert mod.OldV2 is mod.New  # type: ignore[attr-defined]
    finally:
        del sys.modules[mod.__name__]
