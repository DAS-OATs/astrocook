from dataclasses import dataclass, field
from astropy import units as au
import logging
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
import uuid

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SessionMetadataV2:
    """
    Core immutable structure for holding constraints, configuration, and history.

    This structure is designed to be serialized to a JSON or YAML file within
    the .acs archive. It separates the lightweight metadata from the heavy
    spectral arrays.

    Parameters
    ----------
    constraints_by_uuid : dict
        A mapping of Component UUIDs to constraint dictionaries.
    v1_reconstruction_data : Any
        Legacy data structure preserved for technical debt logging or V1 reconstruction.
    log_history_json : Any
        The serialized log object (either a V2 dict or V1 JSON format).
    log_type : str
        Type of log stored: ``'v2'`` (the only one read back), or
        ``'unknown'`` when the session had no log. Older archives may carry
        ``'v1_artifact'`` (a wrapper of V1 logs that was never created,
        removed in phase 6) or ``'v1_legacy'``: their log is not restored.
        Defaults to ``'v1_legacy'``.
    component_uuids : list of str
        An ordered list of all component UUIDs in the session.
    astrocook_version : str or None
        The version of Astrocook that saved the archive; ``None`` for archives
        saved by versions up to 2.0.3, which did not record it (findings
        register F45).
    """
    
    # Stores constraints linked by Component UUID
    constraints_by_uuid: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Stores the original system reconstruction data (for technical debt logging)
    v1_reconstruction_data: Any = None 

    # Improved Session History/Log
    # Stores the serialized log object (either a V2 dict or V1 JSON)
    log_history_json: Any = field(default_factory=dict) # Changed default
    # Stores the type ('v2'; older archives: 'v1_artifact', 'v1_legacy')
    log_type: str = "v1_legacy" # Default to old V1 log
    
    # Stores a list of all component UUIDs in the order they appear
    component_uuids: List[str] = field(default_factory=list)

    # Version of Astrocook that saved the archive (F45)
    astrocook_version: Optional[str] = None

@dataclass(frozen=True)
class DataColumnV2:
    """
    Immutable container for a single column of data with physical units.

    Parameters
    ----------
    values : np.ndarray
        The numerical data array.
    unit : astropy.units.Unit
        The physical unit associated with the values.
    description : str, optional
        A short description of the column's content.
    """
    values: np.ndarray
    unit: au.Unit
    description: str = ""
    
    @property
    def quantity(self) -> au.Quantity:
        """
        Return the column data as an Astropy Quantity.

        Combines the numerical values and the physical unit.
        """
        return au.Quantity(self.values, self.unit)

    # --- [FIX] Compatibility Aliases ---
    @property
    def value(self) -> np.ndarray:
        """
        Return the raw numerical array.

        Alias for ``.values`` to support legacy/Astropy-style access.
        """
        return self.values

    def __len__(self):
        """Return the number of elements in the column."""
        return len(self.values)
    
    def __array__(self):
        """Return the raw array representation (for NumPy compatibility)."""
        return self.values

@dataclass(frozen=True)
class SpectrumDataV2:
    r"""
    Immutable core container for all spectral data.

    This class holds the "heavy" numerical arrays (flux, wavelength, error)
    and auxiliary columns. It is the data payload for :class:`~astrocook.core.spectrum.SpectrumV2`.

    Parameters
    ----------
    x : DataColumnV2
        The independent variable (wavelength or velocity).
    xmin : DataColumnV2
        Lower bound of the x-bins.
    xmax : DataColumnV2
        Upper bound of the x-bins.
    y : DataColumnV2
        The dependent variable (flux density).
    dy : DataColumnV2
        The 1-sigma uncertainty on ``y``.
    aux_cols : dict of DataColumnV2
        Additional columns (e.g. ``'cont'``, ``'model'``, ``'abs_mask'``).
    meta : dict
        Header keywords and other metadata.
    z_em : float
        Emission redshift of the target object.
    resol : float
        Median resolving power :math:`R = \lambda / \Delta\lambda`.
    """
    
    # Dati Spaziali (Assi X e Binning)
    x: DataColumnV2    # Lunghezza d'onda / Velocità
    xmin: DataColumnV2 # Limite inferiore del bin
    xmax: DataColumnV2 # Limite superiore del bin
    
    # Dati Verticali (Flusso e Errore)
    y: DataColumnV2    # Flusso / Densità di flusso
    dy: DataColumnV2   # Errore sul Flusso

    # Colonne Aggiuntive (Cont, Resol, Model, ecc.)
    # Vengono gestite come un dizionario di DataColumnV2
    aux_cols: Dict[str, DataColumnV2] = field(default_factory=dict)
    
    # Metadata di Alto Livello
    meta: Dict[str, Any] = field(default_factory=dict)
    
    # Session properties
    z_em: float = 0.0 # Emission redshift
    resol: float = 0.0 # Median resolution

    @classmethod
    def from_flat_data(cls, 
                       x_values: np.ndarray, x_unit: au.Unit,
                       y_values: np.ndarray, y_unit: au.Unit,
                       dy_values: np.ndarray, 
                       xmin_values: Optional[np.ndarray] = None,
                       xmax_values: Optional[np.ndarray] = None,
                       aux_data: Optional[Dict[str, tuple]] = None, # {'col_name': (values, unit, desc)}
                       meta: Optional[Dict[str, Any]] = None,
                       z_em: float = 0.0,
                       resol: float = 0.0) -> 'SpectrumDataV2':
        """
        Factory method to create a SpectrumDataV2 from raw NumPy arrays.

        Parameters
        ----------
        x_values : np.ndarray
            Wavelength/Velocity data.
        x_unit : astropy.units.Unit
            Unit for the x-axis.
        y_values : np.ndarray
            Flux data.
        y_unit : astropy.units.Unit
            Unit for the y-axis (and dy).
        dy_values : np.ndarray
            Error data.
        xmin_values : np.ndarray, optional
            Lower bin edges. Defaults to ``x_values`` if None.
        xmax_values : np.ndarray, optional
            Upper bin edges. Defaults to ``x_values`` if None.
        aux_data : dict, optional
            Dictionary of auxiliary columns: ``{'name': (values, unit, description)}``.
        meta : dict, optional
            Metadata dictionary.
        z_em : float, optional
            Emission redshift. Defaults to ``0.0``.

        Returns
        -------
        SpectrumDataV2
            The initialized data object.
        """

        # Assicurare che xmin/xmax abbiano valori (es. copiando x se assenti)
        if xmin_values is None: xmin_values = x_values
        if xmax_values is None: xmax_values = x_values

        # Costruire le colonne principali
        x_col = DataColumnV2(x_values, x_unit)
        y_col = DataColumnV2(y_values, y_unit)
        dy_col = DataColumnV2(dy_values, y_unit) # dy usa la stessa unità di y
        xmin_col = DataColumnV2(xmin_values, x_unit)
        xmax_col = DataColumnV2(xmax_values, x_unit)
        
        # Costruire le colonne ausiliarie
        aux_cols = {}
        if aux_data:
            for name, (vals, unit, desc) in aux_data.items():
                 aux_cols[name] = DataColumnV2(vals, unit, desc)
        
        return cls(x=x_col, xmin=xmin_col, xmax=xmax_col, 
                   y=y_col, dy=dy_col, 
                   aux_cols=aux_cols, 
                   meta=meta if meta is not None else {},
                   z_em=z_em,
                   resol=resol)
    
@dataclass(frozen=True)
class ComponentDataV2:
    """
    Immutable data structure for a single absorption component (Voigt profile).
    """
    
    # ID is required
    id: int
    """ Internal integer ID (from V1 legacy logic). """

    # Core Physical Parameters (z, logN, b)
    z: float 
    """ Redshift of the component center """
    dz: Optional[float] 
    """ Uncertainty on redshift. """
    logN: float
    """ Logarithmic column density (cm^-2). """
    dlogN: Optional[float]
    """ Uncertainty on logN. """
    b: float
    """ Doppler broadening parameter (km/s). """
    db: Optional[float]
    """ Uncertainty on b. """
    
    # Turbulance Parameter (from V1)
    btur: float = 0.0
    """ Turbulent broadening contribution (km/s). """
    dbtur: Optional[float] = None
    """ Uncertainty on btur. """

    # Type and identification
    func: str = 'voigt'
    """ Profile function name (default: ``'voigt'``). """
    series: str = 'Ly_a'
    """ Atomic series name (e.g., ``'Ly_a'``, ``'CIV'``). """

    # Fit Metadata
    chi2: Optional[float] = None
    """ Reduced Chi-Squared of the fit group this component belongs to. """
    resol: Optional[float] = None # Resolution used for the fit
    """ Resolution used for the fit. """
    flags: Tuple[str, ...] = ()
    """ Warnings left by the last fit of the component (#517), see
    :data:`astrocook.fitting.voigt_fitter.COMPONENT_FLAGS`: a parameter
    ended at a bound of the fit, or the uncertainties could not be
    computed; or ``not_fitted`` (phase 6-bis): Refit All could not fit its
    group, and it keeps its values from before. Cleared when a parameter is
    edited, and set again by the next fit. """

    # V2 IDENTIFIER: Stable, Global, Immutable
    uuid: str = field(default_factory=lambda: str(uuid.uuid4()))
    """ A globally unique, stable identifier for V2 constraints and linking. """

@dataclass(frozen=True)
class ParameterConstraintV2:
    """
    Immutable structure defining a constraint for a single parameter.

    Used by the fitting engine to determine if a parameter is free, fixed,
    or linked to another component.
    """
    # True if the parameter should be varied during minimization.
    is_free: bool
    """ If True, the parameter varies during minimization. """
    
    # If not None, indicates a link. Must use UUID of the target component.
    target_uuid: Optional[str] = None
    """ If set, this parameter is linked to the component with this UUID. """
    
    # Stores the mathematical expression linking this parameter to the target UUID.
    expression: Optional[str] = None
    """ A mathematical string expression for linking (e.g. ``"p['{uuid}'].z"``). """
    
    # The original V1 integer ID of the target (for debugging/reconstruction only)
    v1_target_id: Optional[int] = None
    """ Legacy ID of the target component (for debugging/reconstruction only). """

@dataclass(frozen=True)
class SystemListDataV2:
    """
    Immutable container for system components and constraints.

    This class replaces the mutable tables used in V1. It holds the list
    of components and the maps defining their relationships.
    """
    
    components: List[ComponentDataV2] = field(default_factory=list)
    """ The list of absorption components. """
    v1_header_constraints: Dict[str, Any] = field(default_factory=dict)
    """ Legacy constraints from V1 file headers. """
    
    # Stores a dictionary mapping (Component UUID, Parameter Name) -> ConstraintDataV2
    parsed_constraints: Dict[Tuple[str, str], ParameterConstraintV2] = field(default_factory=dict)
    """ Map of (Component UUID, Parameter Name) -> ConstraintDataV2. """

    # V2-style constraints: {UUID: {ParamName: {...}}}
    v2_constraints_map: Dict[str, Dict[str, ParameterConstraintV2]] = field(default_factory=dict)
    """ Nested map {UUID: {ParamName: Constraint}}. """
    
    # Add the missing field
    v1_id_to_uuid_map: Dict[int, str] = field(default_factory=dict)
    """ Map of integer IDs to V2 UUIDs. """

    # Placeholder for complex V1 structures (mutable state reference, if needed)
    v1_models_t: Any = None 
    meta: Dict[str, Any] = field(default_factory=dict)
    """ Metadata for the system list. """

@dataclass
class LogEntryV2:
    """
    A single entry in the V2 analysis log.
    
    Corresponds 1-to-1 with a state change in the Session history.
    """
    recipe_name: str
    params: Dict[str, Any]
    #: Why the recipe was run, as given by who ran it (e.g. an agent through
    #: the MCP server); None if not given (phase 6-ter, PR 5a)
    reason: Optional[str] = None


#: The kinds of the notes of a log that :meth:`HistoryLogV2.add_note` takes:
#: ``note`` (what was done or decided outside the recipes), ``gap`` (what the
#: analysis needed and no recipe provided, declared by the agent) and
#: ``gap_seen`` (a gap seen by the MCP server); phase 6-ter, PR 5b
NOTE_KINDS = ("note", "gap", "gap_seen")


@dataclass
class LogNoteV2:
    """
    A note in the V2 analysis log: text that is not a state change (e.g. what
    an agent did outside the recipes), skipped by the replay.

    The notes are kept apart from the entries, so that the entries stay
    1-to-1 with the states of the history, and the logs saved with notes
    replay in the versions before them, which ignore the notes (phase 6-ter,
    PR 5a).
    """
    text: str
    #: The number of entries before the note: it was written at the state
    #: made by the entry ``after - 1`` (0: before any entry)
    after: int
    #: One of :data:`NOTE_KINDS`; a log written by a later version may hold
    #: others, kept as read
    kind: str = "note"

@dataclass
class HistoryLogV2:
    """
    Manages the linear history of operations (Log V2).

    This class maintains a list of :class:`LogEntryV2` objects matching the
    stack of immutable states in the session manager. It supports truncation
    to handle branching history (undo/redo logic).
    """
    # The attributes are documented by the #: comments below and not in an
    # "Attributes" section, which Sphinx would describe twice
    def __init__(self):
        #: The list of operations performed (list of :class:`LogEntryV2`).
        self.entries: List[LogEntryV2] = []
        #: Pointer to the currently active log entry (state), i.e. the last
        #: valid entry.
        self.current_index: int = -1
        #: The notes (list of :class:`LogNoteV2`), in the order they were
        #: written; those written after an entry undone are undone with it.
        self.notes: List[LogNoteV2] = []

    @property
    def is_v2_log(self) -> bool:
        """Return True, indicating this is a V2-native log structure."""
        return True

    def add_entry(self, recipe_name: str, params: Dict[str, Any],
                  reason: Optional[str] = None):
        """
        Add a new log entry, truncating any 'future' history (from Undo actions).

        Parameters
        ----------
        recipe_name : str
            The name of the recipe executed.
        params : dict
            The parameters used for the recipe.
        reason : str, optional
            Why the recipe was run; empty text is no reason.
        """
        # 1. Truncate stale future history, with the notes written in it
        if self.current_index < len(self.entries) - 1:
            logger.debug(f"Truncating V2 log from index {self.current_index + 1}")
            del self.entries[self.current_index + 1:]
        self.notes = [n for n in self.notes if n.after <= self.current_index + 1]
        
        # 2. Add the new entry
        if not (isinstance(reason, str) and reason.strip()):
            reason = None
        entry = LogEntryV2(recipe_name=recipe_name, params=params,
                           reason=reason.strip() if reason else None)
        self.entries.append(entry)
        
        # 3. Update the index to point to the new entry
        self.current_index = len(self.entries) - 1
        logger.debug(f"Added log entry: {recipe_name}. Index is now {self.current_index}.")

    def add_note(self, text: str, kind: str = "note") -> LogNoteV2:
        """
        Add a note at the current entry, without a state.

        The entries undone are kept (a note is not a state, so the redo is
        still possible): the note is undone with the current entry, and
        dropped with it if a new entry replaces it.

        Parameters
        ----------
        text : str
            The text of the note.
        kind : str
            One of :data:`NOTE_KINDS`.

        Returns
        -------
        LogNoteV2
            The note added.

        Raises
        ------
        ValueError
            If the text is empty or the kind is not known.
        """
        if not isinstance(text, str) or not text.strip():
            raise ValueError("A note needs a text.")
        if kind not in NOTE_KINDS:
            raise ValueError(f"Unknown kind of note '{kind}': one of "
                             f"{', '.join(NOTE_KINDS)}.")
        note = LogNoteV2(text=text.strip(), after=self.current_index + 1, kind=kind)
        # After the notes of the same entry or earlier ones, before those of
        # the entries undone
        i = sum(1 for n in self.notes if n.after <= note.after)
        self.notes.insert(i, note)
        return note

    def notes_at(self, after: int) -> List[LogNoteV2]:
        """The notes written after ``after`` entries (at the state made by
        the entry ``after - 1``), in the order they were written."""
        return [n for n in self.notes if n.after == after]

    def timeline(self, include_undone: bool = False) -> List[Tuple[Any, bool]]:
        """
        The entries and the notes in the order they were written.

        Parameters
        ----------
        include_undone : bool
            Also the entries undone and not redone, and their notes.

        Returns
        -------
        list of tuple
            ``(item, active)``: a :class:`LogEntryV2` or a
            :class:`LogNoteV2`, and False if it is undone.
        """
        last = len(self.entries) if include_undone else self.current_index + 1
        items = []
        for i in range(last + 1):
            items += [(n, n.after <= self.current_index + 1) for n in self.notes_at(i)]
            if i < last:
                items.append((self.entries[i], i <= self.current_index))
        return items

    def undo(self) -> bool:
        """
        Move the history pointer back one step.

        Returns
        -------
        bool
            True if undo was successful, False if already before the first
            entry.
        """
        if self.current_index > -1:
            self.current_index -= 1
            logger.debug(f"Log index undone to {self.current_index}")
            return True
        return False

    def redo(self) -> bool:
        """
        Move the history pointer forward one step.

        Returns
        -------
        bool
            True if redo was successful, False if already at the latest entry.
        """
        if self.current_index < len(self.entries) - 1:
            self.current_index += 1
            logger.debug(f"Log index redone to {self.current_index}")
            return True
        return False

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize the active history (entries up to ``current_index``).

        Entries undone and not redone are not part of the saved state and are
        left out.

        Returns
        -------
        dict
            ``{'entries': [{'recipe_name': ..., 'params': ...}, ...],
            'current_index': ...}``; an entry has ``'reason'`` if it was
            given, and the log ``'notes'`` (``[{'text': ..., 'kind': ...,
            'after': ...}, ...]``) if it has any (phase 6-ter, PR 5a): a log
            without either is saved as before.
        """
        active = self.entries[:self.current_index + 1]
        entries = []
        for e in active:
            d = {'recipe_name': e.recipe_name, 'params': e.params}
            if e.reason:
                d['reason'] = e.reason
            entries.append(d)
        out = {'entries': entries, 'current_index': len(active) - 1}
        notes = [{'text': n.text, 'kind': n.kind, 'after': n.after}
                 for n in self.notes if n.after <= len(active)]
        if notes:
            out['notes'] = notes
        return out

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'HistoryLogV2':
        """
        Rebuild a log from the output of :meth:`to_dict`.

        Parameters
        ----------
        data : dict
            The serialized log. Missing or malformed keys give an empty log.

        Returns
        -------
        HistoryLogV2
            The log, pointing at its last entry.
        """
        log = cls()
        for e in (data or {}).get('entries', []) or []:
            log.entries.append(LogEntryV2(recipe_name=e['recipe_name'],
                                          params=dict(e.get('params') or {}),
                                          reason=e.get('reason') or None))
        log.current_index = len(log.entries) - 1
        for n in (data or {}).get('notes', []) or []:
            after = min(max(int(n.get('after', len(log.entries))), 0), len(log.entries))
            log.notes.append(LogNoteV2(text=str(n['text']), after=after,
                                       kind=str(n.get('kind') or 'note')))
        log.notes.sort(key=lambda n: n.after)  # stable: the order of writing
        return log

# --- V1 Compatibility Wrapper ---
