from typing import Any

from astrocook._renamed import renamed
from astrocook.core.structures import SpectrumData, SystemListData

__all__ = ["Session", "SpectrumData", "SystemListData"]

_old_names = renamed(__name__, {"SessionV2": "Session",
                                "SpectrumDataV2": "SpectrumData",
                                "SystemListDataV2": "SystemListData"})


def __getattr__(name: str) -> Any:
    # ``Session`` lives in ``astrocook.session`` since phase 6-quater
    # (PR 5); served lazily so that ``core/`` does not import it.
    if name == "Session":
        from astrocook.session.session import Session
        return Session
    # The old names (PR 6): ``SessionV2`` warns here, then is served as above.
    return _old_names(name)
