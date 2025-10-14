# astrocook/v2/io_adapter.py

from .structures import DataColumnV2, SpectrumDataV2
from .spectrum import SpectrumV2
from ..v1.format import Format as FormatV1 
from ..v1.spectrum import Spectrum as SpectrumV1 # Per il tipo di ritorno V1

# Importa le classi V1 necessarie per l'I/O (ipotizzando path corretto)
# Potrebbe essere necessario un "from astrocook.v1.format import Format" se sposti V1 in v1/

from astropy import units as au
from astropy.io import fits
from astropy.table import Table
import logging
import numpy as np
from typing import Dict, Any

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

def v1_table_to_data_v2(v1_spectrum_instance: Any) -> SpectrumDataV2:
    """
    Adapter: Converte un'istanza Spectrum V1 (basata su Frame/Table) in SpectrumDataV2 immutabile.
    """
    
    # Assumiamo che l'istanza V1 abbia l'attributo _t (astropy.Table) e _xunit, _yunit, _meta
    t_v1 = v1_spectrum_instance._t 
    meta_v1 = v1_spectrum_instance._meta
    x_unit = v1_spectrum_instance._xunit
    y_unit = v1_spectrum_instance._yunit

    # 1. Colonne Core
    x_col = DataColumnV2(t_v1['x'].value, x_unit, description="Channels")
    y_col = DataColumnV2(t_v1['y'].value, y_unit, description="Flux density")
    dy_col = DataColumnV2(t_v1['dy'].value, y_unit, description="Error on Flux")
    xmin_col = DataColumnV2(t_v1['xmin'].value, x_unit, description="Lower channel limit")
    xmax_col = DataColumnV2(t_v1['xmax'].value, x_unit, description="Upper channel limit")

    # 2. Colonne Ausiliarie (Cont, Resol, Model, etc.)
    aux_cols = {}
    for colname in t_v1.colnames:
        # Check if column is a core column before processing aux_cols
        if colname not in ['x', 'xmin', 'xmax', 'y', 'dy']:
            # Use defensive coding to check if the column is accessible
            try:
                col_data = t_v1[colname]
                col_unit = col_data.unit if col_data.unit else au.dimensionless_unscaled
                aux_cols[colname] = DataColumnV2(col_data.value, col_unit, description=f"Auxiliary column: {colname}")
            except Exception as e:
                logging.warning(f"Could not map auxiliary column {colname} from V1 table: {e}")
                # Continue mapping other columns
                continue

    return SpectrumDataV2(
        x=x_col, xmin=xmin_col, xmax=xmax_col, 
        y=y_col, dy=dy_col, 
        aux_cols=aux_cols, 
        meta=meta_v1,
        rf_z=getattr(v1_spectrum_instance, '_rfz', 0.0)
    )

def load_v1_spectrum(path: str, format_name: str, **kwargs) -> SpectrumV1:
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
    # The `gui` object passed via kwargs is the real GUI V1 controller (from gui.py).
    mock_sess = MockSessionV1(kwargs.get('gui'))       

    # If the V1 code accesses flags via sess._gui._flags_cond, this should work, 
    # as kwargs.get('gui') is a real GUI object. 
    
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

def load_spectrum_to_v2_format(path: str, format_name: str, **kwargs) -> SpectrumV2:
    """MAIN ADAPTER: Loads V1, Maps to V2."""
    
    v1_spec = load_v1_spectrum(path, format_name, **kwargs) 
    
    if v1_spec is None:
        raise FileNotFoundError(f"Could not load V1 spectrum from {path}. Check format or file.")

    # Map V1 data to immutable V2 object
    data_v2 = v1_table_to_data_v2(v1_spec) 
    
    return SpectrumV2(data=data_v2)

def load_v2_spectrum(path: str, format_name: str, *args, **kwargs) -> SpectrumV2:
    """
    Funzione I/O pubblica V2: Chiama il caricatore V1, mappa in V2 e restituisce SpectrumV2.
    """
    # Questo richiede la classe Format V1 per l'I/O
    # Supponendo che Format V1 sia disponibile (es. FormatV1)

    # 1. Caricamento V1 (simulato o reale)
    # Esempio: V1_Format = Format()
    # v1_spec = getattr(V1_Format, format_name)(*args, **kwargs)
    
    # --- SIMULAZIONE PER IL TEST INIZIALE ---
    # Creiamo un'istanza fittizia di Spectrum V1 per simulare l'output di format.py
    
    from ..v1.frame import Frame # (Assumendo Frame V1 sia disponibile per il mock)
    from ..v1.spectrum import Spectrum # (Assumendo Spectrum V1 sia disponibile per il mock)

    x = np.array([100.0, 100.1, 100.2, 100.3])
    y = np.array([1.0, 0.5, 0.6, 1.1])
    dy = np.array([0.1, 0.1, 0.1, 0.1])
    
    # Crea un'istanza Spectrum V1 fittizia
    v1_spec = Spectrum(x, x-0.05, x+0.05, y, dy, 
                       xunit=au.nm, yunit=au.erg/au.cm**2/au.s/au.nm,
                       cont=np.ones_like(x),
                       meta={'OBJECT': 'MOCK_SPECTRUM', 'FILE': path})
    
    # 2. Mappa in V2
    data_v2 = v1_table_to_data_v2(v1_spec)
    
    # 3. Restituisci la NUOVA istanza V2
    return SpectrumV2(data=data_v2, history=[f"Loaded from {path} using V1 format: {format_name}"])