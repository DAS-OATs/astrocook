import astropy.units as au
import logging
from typing import TYPE_CHECKING, Optional

from ...v1.vars import xem_d
from ..atomic_data import STANDARD_MULTIPLETS
from ..spectrum import SpectrumV2
from ..structures import HistoryLogV2
try:
    from ...v1.functions import trans_parse
except ImportError:
    logging.error("V1 'trans_parse' not available. Identification recipe will fail.")
    trans_parse = lambda s: []
from ...v1.message import msg_param_fail

if TYPE_CHECKING:
    from ..session import SessionV2

ABSORBERS_RECIPES_SCHEMAS = {
    "identify_lines": {
        "brief": "Identify likely absorption systems.",
        "details": "Finds absorption regions, computes correlation signals for multiplets, and identifies the most likely candidates.",
        "params": [
            {"name": "multiplets", "type": str, "default": "CIV,SiIV,MgII,Ly_a,Ly_b", "doc": "Multiplets to search for (comma-separated)."},
            {"name": "z_grid_dz", "type": float, "default": 1e-4, "doc": "Redshift grid spacing (dz)."},
            {"name": "modul", "type": float, "default": 10.0, "doc": "Significance modulation factor."},
            {"name": "sigma", "type": float, "default": 2.0, "doc": "Minimum S/N peak to be considered."},
            {"name": "distance_pix", "type": int, "default": 3, "doc": "Minimum distance between peaks (pixels)."},
            {"name": "prominence", "type": str, "default": "None", "doc": "Minimum prominence of peaks (e.g., '0.1' or 'None')."},
            {"name": "mask_col", "type": str, "default": "abs_mask", "doc": "Mask column (True=Unabsorbed)"},
            {"name": "min_pix_region", "type": int, "default": 3, "doc": "Minimum width in pixels to count as a region"}
        ],
        "url": "absorbers_cb.html#identify_lines"
    }
}

class RecipeAbsorbersV2:
    def __init__(self, session_v2: 'SessionV2'):
        self._session = session_v2
        self._tag = 'abs'

    def identify_lines(self, 
                         multiplets: str = "CIV,SiIV,MgII,Ly_a,Ly_b",
                         z_grid_dz: str = "1e-4",
                         modul: str = "10.0",
                         sigma: str = "2.0",
                         distance_pix: str = "3",
                         prominence: str = "None",
                         mask_col: str = 'abs_mask', 
                         min_pix_region: str = '3') -> 'SessionV2':
        """
        API: Orchestrator recipe to identify absorption lines.
        (Thin wrapper for SpectrumV2.identify_lines)
        """
        try:
            dz_f = float(z_grid_dz)
            modul_f = float(modul)
            sigma_f = float(sigma)
            dist_i = int(distance_pix)
            prom_f = float(prominence) if prominence.lower() != 'none' else None
            min_pix_i = int(min_pix_region)
            multiplet_list = [m.strip() for m in multiplets.split(',') if m.strip()]
        except ValueError:
            logging.error(msg_param_fail); return 0
        
        try:
            new_spec_v2 = self._session.spec.identify_lines(
                multiplet_list=multiplet_list,
                z_grid_dz=dz_f,
                modul=modul_f,
                sigma=sigma_f,
                distance_pix=dist_i,
                prominence=prom_f,
                mask_col=mask_col,
                min_pix_region=min_pix_i
            )
            return self._session.with_new_spectrum(new_spec_v2)
        except Exception as e:
            logging.error(f"Failed during identify_lines: {e}", exc_info=True)
            return 0