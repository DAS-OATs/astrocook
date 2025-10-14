from .structures import SpectrumDataV2, DataColumnV2

import astropy.units as au
from copy import deepcopy
import numpy as np
from typing import Self, Dict, Any, Optional # Usiamo Self per i metodi che restituiscono una nuova istanza

class SpectrumV2:
    """
    API pubblica e contenitore composito per l'analisi spettrale.
    Tutti i metodi operativi (che trasformano i dati) restituiscono una NUOVA istanza SpectrumV2.
    """
    
    def __init__(self, data: SpectrumDataV2, history: list = None):
        self._data = data
        self.history = history if history is not None else []

    # Accessori (API pubblica: NESSUN setter per i dati)
    @property
    def x(self) -> au.Quantity:
        return self._data.x.quantity

    @property
    def xmin(self) -> au.Quantity:
        return self._data.xmin.quantity

    @property
    def xmax(self) -> au.Quantity:
        return self._data.xmax.quantity

    @property
    def y(self) -> au.Quantity:
        return self._data.y.quantity

    @property
    def dy(self) -> au.Quantity:
        return self._data.dy.quantity

    @property
    def meta(self) -> Dict[str, Any]:
        return deepcopy(self._data.meta) # Restituisce una copia per prevenire modifiche esterne

    @property
    def t(self) -> Dict[str, au.Quantity]:
        """
        Adapter: Restituisce un dizionario simile a un'astropy Table (V1)
        per la compatibilità con la GUI V1 e la visualizzazione.
        """
        # Questa è la funzione ADAPTER per la GUI V1
        output = {
            'x': self.x,
            'xmin': self.xmin,
            'xmax': self.xmax,
            'y': self.y,
            'dy': self.dy,
        }
        # Aggiungere le colonne ausiliarie (es. 'cont', 'resol')
        for name, col_data in self._data.aux_cols.items():
            output[name] = col_data.quantity
            
        return output

    @property
    def _rfz(self) -> float:
        """
        Adapter GETTER per l'attributo V1 '_rfz'.
        Restituisce l'attributo immutabile rf_z dal contenitore dati.
        """
        return self._data.rf_z

    @_rfz.setter
    def _rfz(self, val: float):
        """
        Adapter SETTER per l'attributo V1 '_rfz'.
        Ignora silenziosamente la modifica, poiché SpectrumV2 è immutabile.
        """
        # Se il codice V1 chiama self.spec._rfz = 1.5, l'oggetto V2 non viene modificato.
        # L'unica vera modifica dovrebbe avvenire tramite un metodo che restituisce
        # una NUOVA istanza SpectrumV2 (come shift_rest_frame che avevamo abbozzato).
        pass # Ignora la modifica in-place

    @property
    def _xunit(self) -> au.Unit:
        """
        Adapter GETTER per l'attributo V1 '_xunit'.
        Restituisce l'unità X corrente del contenitore dati.
        """
        return self._data.x.unit
    
    @_xunit.setter
    def _xunit(self, val: au.Unit):
        """
        Adapter SETTER per l'attributo V1 '_xunit'.
        Ignora silenziosamente la modifica (immutabilità V2).
        """
        pass # L'unità X dovrebbe essere aggiornata solo creando un nuovo SpectrumV2
    
    @property
    def _xunit_old(self) -> au.Unit:
        """
        Adapter GETTER per l'attributo V1 '_xunit_old'.
        Questo attributo non esiste nel dato V2, quindi lo emuliamo.
        """
        # Usiamo un fallback nel meta o semplicemente l'unità X stessa.
        # Se il codice V1 ha bisogno di una versione "vecchia" persistente,
        # dovrai modificarne la logica, ma per ora lo emuliamo:
        
        # Tentiamo di recuperare l'unità iniziale dai meta se disponibile, altrimenti l'unità corrente.
        # Se il codice V1 tenta di assegnarlo, il setter vuoto preverrà il crash.
        return self._data.meta.get('initial_x_unit', self._data.x.unit)
    
    @_xunit_old.setter
    def _xunit_old(self, val: au.Unit):
        """
        Adapter SETTER per l'attributo V1 '_xunit_old'.
        Ignora silenziosamente la modifica.
        """
        pass

    @property
    def _yunit(self) -> au.Unit:
        """
        Adapter GETTER per l'attributo V1 '_yunit'.
        Restituisce l'unità Y corrente dal contenitore dati.
        """
        return self._data.y.unit # Restituisce l'unità della colonna y
    
    @_yunit.setter
    def _yunit(self, val: au.Unit):
        """
        Adapter SETTER per l'attributo V1 '_yunit'.
        Ignora silenziosamente la modifica (immutabilità V2).
        """
        pass 

    def _safe(self, quantity: au.Quantity) -> au.Quantity:
        """
        Adapter per il metodo V1 _safe:
        Restituisce una nuova Quantity contenente solo i valori non-NaN.
        (Non modifica lo stato interno come faceva la V1)
        """
        if quantity is None:
            return None
            
        where_safe = ~np.isnan(quantity.value)
        return quantity[where_safe]

    # Metodo generico per accedere a una colonna ausiliaria
    def get_column(self, name: str) -> Optional[au.Quantity]:
        """Accede a una colonna ausiliaria per la GUI (es. 'cont')."""
        col = self._data.aux_cols.get(name)
        return col.quantity if col else None