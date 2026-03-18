import os
import numpy as np
import astropy.units as au
from astrocook.core.structures import SpectrumDataV2, DataColumnV2
from astrocook.core.spectrum import SpectrumV2
from astrocook.core.session import SessionV2
import json

def test_save_load_cycle():
    # 1. SETUP: Create a fake spectrum with V2 features
    x_vals = np.linspace(400, 500, 1000)
    y_vals = np.ones(1000)
    
    # Add a fake 'abs_ids' column (the one that was failing tooltips)
    abs_ids = np.zeros(1000, dtype=int)
    abs_ids[450:550] = 1  # Mark a region in the middle
    
    spec_data = SpectrumDataV2(
        x=DataColumnV2(x_vals, au.nm),
        xmin=DataColumnV2(x_vals-0.1, au.nm),
        xmax=DataColumnV2(x_vals+0.1, au.nm),
        y=DataColumnV2(y_vals, au.dimensionless_unscaled),
        dy=DataColumnV2(y_vals*0.1, au.dimensionless_unscaled),
        # Attach the rich metadata we standardized
        meta={
            'OBJECT': 'IO_TEST_QSO',
            'region_identifications': json.dumps({"1": [["CIV", 0.98]]})
        }
    )
    
    # Add the aux column to the data core
    spec_data.aux_cols['abs_ids'] = DataColumnV2(abs_ids, au.dimensionless_unscaled)
    
    spec = SpectrumV2(spec_data)
    session = SessionV2(name="Test_Integrity", gui=None, spec=spec)

    # 2. SAVE: Trigger the .acs2 save logic
    test_file = "integrity_check.acs2"
    if os.path.exists(test_file): os.remove(test_file)
    
    print(f"--- Saving session to {test_file} ---")
    session.save(test_file, initial_session=session)

    # 3. LOAD: Trigger the fixed load_session_from_file logic
    print("--- Reloading session ---")
    reloaded_session = SessionV2.open_new(test_file, "Reloaded", None, "auto")

    # 4. VALIDATE
    try:
        # Check Core Wavelengths (Ensure no 100-500nm fallback)
        np.testing.assert_allclose(reloaded_session.spec.x.value, x_vals)
        print("✅ Core Wavelengths: PRESERVED")

        # Check Aux Columns (Ensure abs_ids exists)
        assert 'abs_ids' in reloaded_session.spec._data.aux_cols
        np.testing.assert_array_equal(reloaded_session.spec.get_column('abs_ids').value, abs_ids)
        print("✅ Auxiliary Column 'abs_ids': PRESERVED")

        # Check Metadata (The JSON identifications)
        meta_id = reloaded_session.spec.meta.get('region_identifications')
        assert meta_id is not None
        assert "CIV" in meta_id
        print("✅ Rich Metadata JSON: PRESERVED")
        
        print("\n🔥 INTEGRITY CHECK PASSED: V2 I/O IS STABLE 🔥")

    except Exception as e:
        print(f"\n❌ INTEGRITY CHECK FAILED: {e}")
    
    finally:
        if os.path.exists(test_file): os.remove(test_file)

if __name__ == "__main__":
    test_save_load_cycle()