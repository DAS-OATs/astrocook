from typing import List

from .structures import ComponentDataV2, SystemListDataV2

# --- 3. System List API Layer (The orchestrator) ---
class SystemListV2:
    """API layer for managing component lists, constraints, and fitting."""
    
    def __init__(self, data: SystemListDataV2):
        self._data = data
        # NOTE: ConstraintModelV2 will be instantiated here later
        # self.constraints = ConstraintModelV2(data) 
        
    @property
    def components(self) -> List[ComponentDataV2]:
        return self._data.components
        
    # TODO: Methods for constraints, fitting, and immutable updates will go here.