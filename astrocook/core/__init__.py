from typing import Any

from astrocook.core.structures import SpectrumDataV2, SystemListDataV2

__all__ = ["SessionV2", "SpectrumDataV2", "SystemListDataV2"]


def __getattr__(name: str) -> Any:
    # ``SessionV2`` lives in ``astrocook.session`` since phase 6-quater
    # (PR 5); served lazily so that ``core/`` does not import it.
    if name == "SessionV2":
        from astrocook.session.session import SessionV2
        return SessionV2
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
