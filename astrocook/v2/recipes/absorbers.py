import logging
from typing import TYPE_CHECKING, Optional
from ..structures import HistoryLogV2
from ...v1.message import msg_param_fail

if TYPE_CHECKING:
    from ..session import SessionV2

ABSORBERS_RECIPES_SCHEMAS = {
    "detect_regions": {
        "brief": "Detect absorption regions.",
        "details": "Group consecutive absorbed pixels (from 'mask_unabs') into labeled regions.",
        "params": [
            {"name": "mask_col", "type": str, "default": "mask_unabs", "doc": "Continuum mask column (True=Unabsorbed)"},
            {"name": "min_pix", "type": int, "default": 3, "doc": "Minimum width in pixels to count as a region"}
        ],
        "url": "absorbers_cb.html#detect_regions"
    }
}

class RecipeAbsorbersV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'abs'

    def detect_regions(self, mask_col: str = 'mask_unabs', min_pix: str = '3') -> 'SessionV2':
        """
        API: Detects absorption regions and creates 'region_id' column.
        """
        try:
            min_pix_i = int(min_pix)
        except ValueError:
            logging.error(msg_param_fail); return 0
            
        try:
            new_spec_v2 = self._session.spec.detect_absorptions(
                mask_col=mask_col,
                min_pix=min_pix_i
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during detect_regions: {e}", exc_info=True)
            return 0