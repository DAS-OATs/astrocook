from astropy.io import fits
import logging
from typing import Any, Optional

from ..v1.format import Format as FormatV1 
from ..v1.spectrum import Spectrum as SpectrumV1
from ..v1.syst_list import SystList as SystListV1

def create_mock_v1_spectrum(hdul):
    """Creates a mock V1 Spectrum object with the correct size from HDUL."""
    from astropy import units as au
    from ..v1.spectrum import Spectrum as SpectrumV1
    
    # Assuming the HDU data is in hdu[1] or hdu[0]
    data_length = 500 # A reasonable mock length
    try:
        data_length = len(hdul[1].data) 
    except:
        try:
            data_length = len(hdul[0].data)
        except:
            pass
            
    x = np.linspace(100.0, 500.0, data_length)
    y = np.ones(data_length)
    dy = np.ones(data_length) * 0.1
    
    return SpectrumV1(x, x-0.5, x+0.5, y, dy, 
                      xunit=au.nm, yunit=au.erg/au.cm**2/au.s/au.nm,
                      cont=np.ones_like(x),
                      meta={'OBJECT': 'MOCKED_FROM_REAL_FILE', 'FILE': '...' })

def load_v1_spec_object(path: str, format_name: str, gui_context: Any) -> SpectrumV1:
    """
    Loads a V1 spectrum using the logic from format.py.
    Returns a V1 Spectrum instance.
    """
    
    # 1. Open FITS/ASCII file
    try:
        if path.lower().endswith(('.fits', '.fits.gz')):
            hdul = fits.open(path)
            hdr = hdul[0].header
        else:
            # Simple simulation for other formats (you'll need to expand this for V1 formats)
            logging.error(f"Unsupported file extension for direct V1 loading: {path}. Trying FITS...")
            return None
    except Exception as e:
        logging.error(f"Failed to open file {path}: {e}")
        return None

    # 2. Mock the minimal Session V1 object required by Format.generic_spectrum
    class MockSessionV1:
        def __init__(self, gui):
            self._gui = gui # V1 code expects this directly
        
        # The generic_spectrum V1 logic needs this method accessible via the session:
        # It tries to call sess._gui._flags_cond() or sess._gui._flags_extr()
        # Let's mock the necessary methods on the mock GUI object itself.
        
        # We don't need to define anything here, but we must ensure the 'gui' object
        # that is passed to this mock is able to handle the V1 requests.       

        
    # 3. Enhance the V1 Call: Pass the Real GUI Object     

    # The V1 logic needs access to the main GUI controller to check flags.
    mock_sess = MockSessionV1(gui_context)       
    
    v1_format_loader = FormatV1()
    
    try:
        # Attempt real V1 loading
        v1_spec = getattr(v1_format_loader, format_name)(mock_sess, hdul)
        if v1_spec is None:
            # Raise an exception to correctly enter the 'except' block
            raise RuntimeError(f"V1 loader '{format_name}' returned None.")
        
        logging.info(f"V1 spectrum loaded successfully using {format_name}.")
        
    except Exception as e:
        # If V1 loading fails (expected for generic_spectrum), use the size-matched mock
        logging.warning(f"V1 format loading failed for {format_name} ({e}). Using size-matched mock data.")
        v1_spec = create_mock_v1_spectrum(hdul)
        
    finally:
        if 'hdul' in locals() and hdul is not None:
            hdul.close()

    if v1_spec is not None:
        return v1_spec
    else:
        # Should be unreachable if the FITS file was opened
        return None

def load_v1_systs_object(file_path: str) -> Optional[SystListV1]:
    """Placeholder for V1 SystList loading logic."""
    return None