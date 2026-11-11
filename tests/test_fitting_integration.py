import pytest
import numpy as np
from astropy.table import Table
from astrocook.fitting.voigt_model import VoigtModelConstraint
from astrocook.core.structures import SystemListData
from astrocook.core.system_list import SystemList
from astrocook.core.structures import ComponentData

# --- Helper function to simulate a V2 System List with V1 Model Placeholder ---
def create_v2_system_list_with_v1_models():
    """
    Creates a SystemList object containing one ComponentData and a mock V1 _mods_t table.
    """
    # 1. Create a minimal ComponentData
    comp_data = ComponentData(id=1, z=2.5, dz=0.0, logN=14.0, dlogN=0.1, b=10.0, db=1.0)
    
    # 2. Create a mock V1 _mods_t structure (an Astropy Table)
    mock_v1_models_t = Table({'id': [[1], [2]], 'mod': [object(), object()]})
    
    # 3. Create the SystemListData core
    data_core = SystemListData(
        components=[comp_data],
        v1_models_t=mock_v1_models_t, # The placeholder for the lmfit models
        meta={}
    )
    
    # 4. Create and return the API layer object
    return SystemList(data=data_core)

# --- The Integration Test ---
def test_constraint_model_loads_v1_placeholder():
    system_list_v2 = create_v2_system_list_with_v1_models()
    
    # 1. Instantiate the Constraint Model
    # The constructor should successfully pull the v1_models_t from the system list.
    try:
        constraint_model = VoigtModelConstraint(system_list_v2)
        
        # 2. Assert successful integration of the placeholder data
        # We check that the model object has stored the placeholder data internally.
        assert constraint_model._system_list is system_list_v2
        
        # 3. Assert the placeholder data is accessible (for future constraint parsing)
        # We access the V1 placeholder via the system_list's data core
        placeholder_ref = constraint_model._system_list.v1_models_t
        assert isinstance(placeholder_ref, Table)
        assert len(placeholder_ref) == 2
        
    except Exception as e:
        pytest.fail(f"VoigtModelConstraint initialization failed: {e}")

def test_p_full_vector_initialization():
    system_list_v2 = create_v2_system_list_with_v1_models()
    
    # 1. Instantiate the Constraint Model
    constraint_model = VoigtModelConstraint(system_list_v2)
    
    # 2. Retrieve the P_full vector
    p_full_vector = constraint_model._initial_p_vector
    
    # 3. Define the expected full vector (z, logN, b, btur for one component)
    # Based on MockSystV1: z=2.5, logN=14.0, b=10.0, btur=0.0
    expected_full = np.array([2.500000, 14.0, 10.0, 0.0])
    
    # 4. Assert length and values
    assert len(p_full_vector) == 4 # We defined 4 parameters in PARAMETER_ORDER
    
    # Assert that the initial vector matches the expected values
    assert np.allclose(p_full_vector, expected_full)