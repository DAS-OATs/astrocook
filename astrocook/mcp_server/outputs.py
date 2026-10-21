"""Structured outputs of the tools.

The SDK derives the output schema of each tool from these models, and the
client shows it to the model next to the tool description: the field
descriptions are part of the interface.
"""
from typing import Optional

import numpy as np
from pydantic import BaseModel, Field

from astrocook.core.session import SessionV2


class SpectrumSummary(BaseModel):
    """What the agent needs to know about a spectrum, without the data."""

    n_pixels: int = Field(description="Number of pixels.")
    x_min: float = Field(description="First pixel centre, in x_unit.")
    x_max: float = Field(description="Last pixel centre, in x_unit.")
    x_unit: str = Field(description="Unit of the x axis (e.g. nm).")
    y_unit: str = Field(description="Unit of the flux.")
    z_em: float = Field(description="Emission redshift; 0 if not set.")
    resol: float = Field(description=(
        "Resolving power R stored with the spectrum; 0 if not set."))
    has_continuum: bool = Field(description="Whether a continuum is present.")
    n_components: int = Field(description="Number of absorption components.")


class SessionInfo(BaseModel):
    """A session after opening it or running a recipe on it."""

    session_id: str = Field(description=(
        "Identifier of the session, to pass to the other tools."))
    spectrum: SpectrumSummary
    log_entries: int = Field(description=(
        "Number of operations recorded in the session log."))
    saved_with_version: Optional[str] = Field(None, description=(
        "Astrocook version that saved the archive, if opened from one "
        "that records it."))
    warnings: list[str] = Field(default_factory=list, description=(
        "Warnings logged by Astrocook during the call."))


class SaveResult(BaseModel):
    """Where a session was saved."""

    session_id: str
    path: str = Field(description="Absolute path of the .acs2 archive.")
    log_entries: int = Field(description=(
        "Number of operations recorded in the saved log."))
    warnings: list[str] = Field(default_factory=list, description=(
        "Warnings logged by Astrocook during the call."))


def spectrum_summary(session: SessionV2) -> SpectrumSummary:
    """Summarize the spectrum and the components of a session."""
    spec = session.spec
    x = spec.x
    return SpectrumSummary(
        n_pixels=len(x),
        x_min=float(np.nanmin(x.value)),
        x_max=float(np.nanmax(x.value)),
        x_unit=x.unit.to_string(),
        y_unit=spec.y.unit.to_string(),
        z_em=float(spec._data.z_em),
        resol=float(spec._data.resol),
        has_continuum=spec.cont is not None,
        n_components=len(session.systs.components) if session.systs else 0,
    )
