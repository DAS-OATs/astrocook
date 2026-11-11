"""Typing helper of the recipes (phase 6-quater, PR 3)."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from astrocook.session.session import SessionV2
    from astrocook.core.spectrum import SpectrumV2


def spec_of(session: SessionV2) -> SpectrumV2:
    """The spectrum of a session, for the type checker.

    ``SessionV2.spec`` is typed ``Optional`` (a session may hold no
    spectrum), but a recipe that works on the spectrum runs on a session
    that has one: this narrows the type, with no change of behaviour.
    """
    spec = session.spec
    assert spec is not None, "the session has no spectrum"
    return spec


def column_of(spec: SpectrumV2, name: str) -> Any:
    """A column that the caller knows is there (``get_column`` is typed to
    return None for a missing one)."""
    col = spec.get_column(name)
    assert col is not None, f"the spectrum has no column {name!r}"
    return col
