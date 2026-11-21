from astrocook._renamed import renamed
import logging
import os
import numpy as np
from astropy.io import ascii
import astropy.units as au
from astropy.table import Column, Table
from typing import TYPE_CHECKING, List, Union

from astrocook.core.column_names import internal_column, internal_expression, is_axis_column
from astrocook.core.units import known_unit
from astrocook.recipes.errors import InvalidParameter, RecipeFailed
from astrocook.recipes._spectrum import spec_of
from astrocook.recipes.params import typed_recipe

if TYPE_CHECKING:
    from astrocook.session.session import Session

logger = logging.getLogger(__name__)


# --- The text of the methods (phase 6-ter, PR 5c: astrocook.methods_text) ---
from astrocook.methods_text import file_name
def _methods_import_ascii(p):
    return (f"Columns were read from the table {file_name(p, 'file_path')} and added "
            f"to the spectrum{p.menu}")


def _methods_import_model(p):
    source = (file_name(p, 'file_path') if p.raw('file_path')
              else f"the session {p.raw('source_session')}")
    column = f" (column {p.raw('column')})" if p.raw('column') else ""
    words = {'telluric': "telluric transmission",
             'cont': "continuum"}.get(str(p.raw('role')), "model")
    return (f"An external {words} was imported from {source}{column} and brought "
            f"onto the grid of the spectrum (copied on the same grid, otherwise "
            f"averaged over its pixels){p.menu}")


FILE_RECIPES_SCHEMAS = {
    "export_ascii": {
        "methods": "The data ({targets}) were exported as CSV files{menu}",
        "brief": "Export session data to ASCII.",
        "url": "tutorials/getting_started.html#export-to-ascii",
        "menu": "File > &Export to ASCII",
        "details": "Writes the chosen data as CSV files in a folder: the columns of the spectrum to <prefix>_spec.csv (always with the wavelength), the system list to <prefix>_systems.csv (always with the series and the redshift). Items that are neither a column nor 'spec' or 'systems' are skipped, with a warning.",
        "params": [
            {"name": "targets", "type": str, "default": "F, cont, systems", "doc": "Comma-separated items to export: 'spec' (all the columns of the spectrum) or some of its columns (e.g. 'F, dF, cont'; the file names them as in the data: x, y, dy), and 'systems' (the system list) or some of its columns (e.g. 'z, logN')."},
            # 'path' until phase 6-quater, PR 7 (an agent could take it for
            # the path of a file, as in session_save)
            {"name": "folder", "type": str, "default": "", "gui_hidden": True, "aliases": ["path"], "doc": "Folder where the files are written."},
            # F106: given by the file dialog of the GUI, required elsewhere
            {"name": "prefix", "type": str, "default": "", "gui_hidden": True, "doc": "Beginning of the file names: <prefix>_spec.csv and <prefix>_systems.csv."},
            # F106: the GUI asks before and passes True; a script or an
            # agent says it (the files were overwritten without a word)
            {"name": "overwrite", "type": bool, "default": False, "gui_hidden": True, "doc": "Replace the files if they exist; if False, the export stops and names them."}
        ],
    },
    "import_ascii": {
        "methods": _methods_import_ascii,
        "brief": "Import session data from ASCII.",
        "url": "tutorials/getting_started.html#import-from-ascii",
        "menu": "File > Import from &ASCII",
        "details": "Reads columns from an ASCII/CSV table with a header and adds them to the spectrum, replacing those with the same name; the table must have one row per pixel of the spectrum. A column x is not read: the grid of the spectrum is kept.",
        "params": [
            {"name": "file_path", "type": str, "default": "", "gui_hidden": True, "doc": "Path of the ASCII/CSV file to read."}
        ],
    },
    "import_model": {
        "methods": _methods_import_model,
        "brief": "Import an external model into the spectrum.",
        "url": "tutorials/external_models.html#importing-a-model",
        "menu": "File > Import &Model",
        "details": "Reads a model (a telluric transmission, a continuum) from a file or another session and writes it in the column of its role, on the grid of the spectrum: copied on the same grid, otherwise averaged over the pixels. A telluric model is 1 beyond its range, a continuum NaN.",
        "params": [
            {"name": "file_path", "type": str, "default": "", "doc": "Path of the file of the model: a text table with named columns, the CSV of Export to ASCII, a FITS table. Empty if the model is in another session."},
            {"name": "source_session", "type": str, "default": "", "kind": "session", "doc": "Name of the open session that holds the model (in a script, also the session itself); empty if it is in a file."},
            # T4 adds 'model' (recipes/_models.py)
            {"name": "role", "type": str, "default": "telluric", "choices": ["telluric", "cont"], "doc": "Role of the model: 'telluric', a transmission written in the column telluric; 'cont', a continuum in the unit of the flux written in the column cont (cont × telluric is the normalization)."},
            {"name": "column", "type": str, "default": "", "doc": "Column of the model in the source; empty: found by name (telluric: telluric, telluric_model, trans_molecfit, mtrans; cont: cont, continuum)."},
            # T3 (Guido, 20 November 2026): as in fit_continuum
            {"name": "renorm_model", "type": bool, "default": True, "doc": "With role 'cont', rescale the model to the new continuum (asked by the GUI).", "gui_hidden": True, "agents": False}
        ],
    }
}

#: The names of the columns of every spectrum in the CSV files: those of the
#: GUI, in ASCII (phase 6-quater, PR 8-bis; x, xmin, xmax, y, dy before)
_CSV_NAMES = {'x': 'lam', 'xmin': 'lammin', 'xmax': 'lammax', 'y': 'F', 'dy': 'dF'}


def _spectrum_table(spec, columns) -> Table:
    """The columns of the spectrum as a table for a CSV file: the names of
    the GUI in ASCII, the units, and z_em in the metadata (F129: the files
    had neither, and the loader assumed nm and erg/(Angstrom s cm2))."""
    table = Table()
    for name in columns:
        q = spec.t[name]
        unit = getattr(q, 'unit', None)
        table[_CSV_NAMES.get(name, name)] = Column(
            np.asarray(getattr(q, 'value', q)),
            unit=None if unit == au.dimensionless_unscaled else unit)
    table.meta['z_em'] = float(spec._data.z_em)
    return table


class RecipeFile:
    """
    Recipes for file input/output operations.

    This class manages exporting session data (spectra and system lists) to 
    external formats like ASCII/CSV, and importing modified data back into 
    the session. It operates immutably where applicable to support Undo/Redo.
    """
    def __init__(self, session_v2: 'Session'):
        self._session = session_v2
        self._tag = 'cb'

    @typed_recipe
    def export_ascii(self, targets: str, folder: str, prefix: str,
                     overwrite: bool = False) -> bool:
        """
        Export selected components of the session to ASCII/CSV files.

        Generates distinct files for spectrum data (`_spec.csv`) and system 
        lists (`_systems.csv`) based on the requested targets. Spectrum 
        exports automatically include the wavelength array (``lam``), while 
        system exports automatically include the transition name (`series`) 
        and redshift (`z`).

        The spectrum is written as ECSV (the CSV format of astropy, phase
        6-quater, PR 8-bis, F129): a header of comment lines with the unit
        of each column and the emission redshift, then the names and the
        values separated by commas (``pandas.read_csv(path, comment='#')``
        reads it). The columns have the names of the GUI in ASCII: ``lam``,
        ``lammin``, ``lammax``, ``F``, ``dF`` and the others by their own
        name (before: ``x``, ``xmin``, ``xmax``, ``y``, ``dy``, without
        units). Opened again, the file gives the same columns and units.

        Parameters
        ----------
        targets : str
            Comma-separated list of structures or columns to export 
            (e.g., ``'spec, systems'`` or ``'F, cont, logN'``).
        folder : str
            The folder where the files are written (``path`` until phase
            6-quater, PR 7: the old name is still accepted).
        prefix : str
            The base filename prefix for the exported files.
        overwrite : bool, optional
            Replace the files if they exist. If ``False`` (the default) and
            a file exists, nothing is written and the recipe fails naming the
            files (F106: they were overwritten without a word outside the
            GUI, which asks before and passes ``True``).

        Returns
        -------
        bool
            ``True`` if at least one file was saved, ``False`` if the targets
            name nothing to export.

        Raises
        ------
        RecipeError
            If a file cannot be written (phase 6-bis, PR 3: it returned False).
        """
        spec_req, write_systems, syst_req, target_list = self._export_targets(targets)
        spec_colnames = spec_of(self._session).t.colnames if self._session.spec else []
        syst_colnames = self._session.systs.t.colnames if self._session.systs else []
        skipped = [t for t in target_list if t not in ('spec', 'systems')
                   and t not in spec_colnames and t not in syst_colnames]
        if skipped:
            # Phase 6-quater, PR 7: they were skipped without a word (e.g. a
            # misspelt 'systems', or the default 'cont' with no continuum)
            logger.warning(f"Not exported, not found in the session: {', '.join(skipped)}.")
        outputs = (['spec'] if spec_req else []) + (['systems'] if write_systems else [])
        existing = [f"{prefix}_{what}.csv" for what in outputs
                    if os.path.exists(os.path.join(folder, f"{prefix}_{what}.csv"))]
        if existing and not overwrite:
            raise RecipeFailed(f"{', '.join(existing)} already "
                               f"{'exists' if len(existing) == 1 else 'exist'} in {folder}: "
                               f"nothing is written. Use overwrite=True to replace "
                               f"{'it' if len(existing) == 1 else 'them'}.")
        saved_any = False

        if spec_req:
            spec = spec_of(self._session)
            out_path = os.path.join(folder, f"{prefix}_spec.csv")
            ascii.write(_spectrum_table(spec, ['x'] + spec_req), out_path,
                        format='ecsv', delimiter=',', overwrite=True)
            saved_any = True

        # 2. Systems Export (Refined)
        if write_systems:
            syst_table = self._session.systs.t
            
            if "systems" not in target_list:
                # [NEW LOGIC] Always include 'series' and 'z' for context.
                # We use dict.fromkeys to maintain order and uniqueness.
                # 'series' first, then 'z', then other requested params.
                base_cols = ['series', 'z']
                # Add any other columns requested that aren't already z/series
                extra_cols = [t for t in syst_req if t not in base_cols]
                
                cols_to_keep = base_cols + extra_cols
                # Verify they exist in the actual table before slicing
                cols_to_keep = [c for c in cols_to_keep if c in syst_table.colnames]
                syst_table = syst_table[cols_to_keep]
            
            out_path = os.path.join(folder, f"{prefix}_systems.csv")
            ascii.write(syst_table, out_path, format='csv', overwrite=True)
            saved_any = True
            
        return saved_any

        
    def export_ascii_outputs(self, targets: str) -> List[str]:
        """
        The files that :meth:`export_ascii` writes for these targets.

        The GUI asks this before exporting, to suggest a file name and to
        check for files that would be overwritten (F82: it used to decide
        it by itself, without the names of the columns shown in the GUI).

        Parameters
        ----------
        targets : str
            Comma-separated list of structures or columns, as for
            :meth:`export_ascii`.

        Returns
        -------
        list of str
            ``'spec'`` (the ``_spec.csv`` file) and/or ``'systems'`` (the
            ``_systems.csv`` file), in this order.
        """
        spec_req, write_systems, _, _ = self._export_targets(targets)
        return (['spec'] if spec_req else []) + (['systems'] if write_systems else [])

    def _export_targets(self, targets: str):
        """What :meth:`export_ascii` writes for these targets.

        Returns the columns of the spectrum to write (besides ``x``, always
        written), whether the systems are written, the columns of the systems
        asked for, and the targets with the names of the GUI translated.
        ``'spec'`` asks for all the columns of the spectrum.
        """
        # The names of the GUI (λ, F, dF, ...) to those of the data (F120)
        processed_targets = internal_expression(targets)

        target_list = [t.strip() for t in processed_targets.split(',') if t.strip()]

        # Get authoritative column lists
        spec_colnames = spec_of(self._session).t.colnames if self._session.spec else []
        syst_colnames = self._session.systs.t.colnames if self._session.systs else []

        # Any target that exists in spec table (and isn't the literal 'systems')
        if 'spec' in target_list:
            spec_req = [c for c in spec_colnames if c != 'x']
        else:
            spec_req = [t for t in target_list if t in spec_colnames and t != 'x']
        syst_req = [t for t in target_list if t in syst_colnames]
        write_systems = bool(("systems" in target_list or syst_req) and self._session.systs)
        return spec_req, write_systems, syst_req, target_list

    @typed_recipe
    def import_ascii(self, file_path: str) -> 'Session':
        """
        Import columns from an ASCII/CSV file into the current session.

        Reads external data and overwrites or appends columns in the 
        spectrum based on matching column headers. The imported file must 
        have the exact same number of rows as the current spectrum to ensure 
        grid alignment. 

        Parameters
        ----------
        file_path : str
            The full path to the ASCII/CSV file to import.

        Returns
        -------
        Session
            A new session instance containing the imported data.

        Raises
        ------
        RecipeError
            If the file cannot be read, or its length is not that of the
            spectrum.
        """
        table = ascii.read(file_path)
        current_len = len(spec_of(self._session).x)
        import_len = len(table)

        # 2. Safety Check: Pixel count must match
        if import_len != current_len:
            raise RecipeFailed(f"The file has {import_len} pixels and the session "
                               f"{current_len}: they must be the same.")

        new_spec = spec_of(self._session)
        self._check_ascii_axis(table, new_spec, file_path)
        imported_cols = []
        replaced = []

        # 3. Iterate through columns in the ASCII file
        for file_name_col in table.colnames:
            # The names of the GUI (F, dF, lam, ...) are those of the data
            # (PR 8-bis: F made a column F)
            col_name = internal_column(file_name_col)
            # The axis is skipped to prevent grid misalignment (x, and since
            # PR 8-bis lam, xmin, xmax and their other names)
            if col_name.lower() == 'x' or is_axis_column(col_name):
                continue

            values = np.array(table[file_name_col])
            unit = table[file_name_col].unit
            if col_name in ('y', 'dy') and unit is not None:
                unit = known_unit(unit)
                if unit != new_spec.y.unit:
                    raise RecipeFailed(f"The column {file_name_col} of the file is in {unit}, "
                                       f"the flux of the session in {new_spec.y.unit}: "
                                       "they must be the same.")
            # Update the column (this handles both core y/dy and aux
            # columns); the unit of the file, if it has one (PR 8-bis)
            if col_name in new_spec.t.colnames:
                replaced.append(file_name_col)
            new_spec = new_spec.update_column(col_name, values,
                                              unit=None if unit is None else known_unit(unit))
            imported_cols.append(file_name_col)

        logger.info(f"Columns imported: {', '.join(imported_cols)}.")
        if replaced:
            # FT3 (Guido, 20 November 2026): they were replaced without a word
            logger.warning(f"Columns of the spectrum replaced by those of the file: "
                           f"{', '.join(replaced)}.")
        return self._session.with_new_spectrum(new_spec)

    @staticmethod
    def _check_ascii_axis(table, spec, file_path: str) -> None:
        """FT3 (Guido, 20 November 2026): a table whose axis is not that of
        the spectrum is refused (only the number of rows was checked). The
        axis of the file (``x``, ``lam``, ...; in nm without a unit) must be
        within a tenth of a pixel of the axis of the spectrum, as for the
        copy of ``import_model``; a table without an axis is taken as it
        is."""
        from astrocook.core.model_grid import same_grid
        from astrocook.io.tables import AXIS_NAMES, find_column
        name = find_column(table, AXIS_NAMES)
        if name is None:
            return
        unit = known_unit(table[name].unit)
        if unit is None or unit == au.dimensionless_unscaled:
            unit = au.nm
        try:
            x_file = (np.asarray(table[name], dtype=float) * unit).to(
                spec.x.unit, equivalencies=au.spectral()).value
        except (au.UnitConversionError, TypeError, ValueError):
            x_file = None
        if x_file is None or not same_grid(x_file, spec.x.value, spec.xmin.value,
                                           spec.xmax.value):
            raise RecipeFailed(
                f"The axis of {os.path.basename(file_path)} ({name}) is not that of the "
                "spectrum: the columns would go to other pixels. To bring a model "
                "onto the grid of the spectrum, use Import Model (import_model).")

    @typed_recipe
    def import_model(self, file_path: str = "", source_session: Union[str, 'Session'] = "",
                     role: str = "telluric", column: str = "",
                     renorm_model: bool = True) -> 'Session':
        """
        Import an external model into the spectrum (phase 7-ter, PR T2).

        A model is a column of the spectrum with a role. The model is read
        from a file or from another session, brought onto the grid of the
        spectrum and written in the column of its role, replacing it if it
        exists (undo, log and archive as for any recipe):

        - the grid: a copy when the axis of the model is that of the
          spectrum within a tenth of a pixel, otherwise the average of the
          model over each pixel of the spectrum (``xmin``-``xmax``), each
          pixel of the model weighted by its overlap
          (:mod:`astrocook.core.model_grid`); a pixel the model does not
          cover in full gets the value of the role beyond the model;
        - the role ``telluric``: a transmission (dimensionless, 0-1)
          written in the column ``telluric``, 1 beyond the model; the
          normalization is ``cont × telluric`` unless the flux is flagged
          as corrected (``set_properties(telluric_corrected=...)``);
        - the role ``cont`` (T3): a continuum, in the unit of the flux,
          written in the column ``cont``, NaN beyond the model (there the
          flux is not normalized; a WARNING says on how many pixels). A
          unit of the source is converted to that of the flux, refused if
          it cannot be; without a unit the values are taken in the unit of
          the flux, and refused when their median and that of the flux,
          over the pixels covered, differ by more than a factor 10 (a
          normalized continuum under a physical flux). The ``model`` of
          the spectrum follows the new continuum if ``renorm_model``, as in
          ``fit_continuum`` (kept as it was where there is no continuum).

        The models are taken in the frame of the spectrum: no shift. A file
        is read by the names of its columns (:mod:`astrocook.io.tables`):
        a text table separated by spaces, tabs or commas with a header,
        also commented (``# x xmin xmax trans_molecfit dy``), the CSV of
        ``export_ascii``, a FITS table (Astrocook, the Molecfit outputs
        with ``mtrans``); the axis is ``x``, ``lam``, ``lambda``, ``wave``
        or ``wavelength`` (in nm in a text file without units), its edges
        ``xmin``-``xmax`` (``lammin``-``lammax``) if the file has them. The
        other columns are ignored (the ``dy`` of a Molecfit model is the
        error of the science spectrum, not of the model).

        Parameters
        ----------
        file_path : str
            The file of the model; empty if it is in another session.
        source_session : str or Session
            The open session that holds the model (its name; in a script
            also the session itself); empty if it is in a file.
        role : str
            ``'telluric'`` or ``'cont'`` (the role ``'model'`` comes in
            phase 7-ter, T4).
        column : str
            The column of the model in the source; empty: found by the
            names of the role (for ``'telluric'``: ``telluric``,
            ``telluric_model``, ``trans_molecfit``, ``mtrans``; for
            ``'cont'``: ``cont``, ``continuum``; ignoring the case).
        renorm_model : bool
            With the role ``'cont'``, rescale the ``model`` of the spectrum
            to the new continuum (the GUI asks); ignored by the other roles.

        Returns
        -------
        Session
            A new session with the model in its column.

        Raises
        ------
        RecipeError
            If neither or both sources are given, the source cannot be read
            or has no column of the model, its values cannot be this model,
            or the model does not overlap the spectrum.
        """
        from astrocook.recipes import _models
        from astrocook.recipes._sessions import open_session
        has_file = bool(str(file_path).strip())
        has_session = bool(source_session) and (not isinstance(source_session, str)
                                                or bool(source_session.strip()))
        if has_file == has_session:
            raise InvalidParameter("import_model: give file_path or source_session "
                                   "(one of the two).")
        model_role = _models.MODEL_ROLES[role]
        spec = spec_of(self._session)
        if has_file:
            source = _models.source_from_file(str(file_path).strip(), model_role,
                                              column.strip())
        else:
            other = open_session(self._session, source_session, 'source_session')
            source = _models.source_from_session(other, model_role, column.strip())
        return self._session.with_new_spectrum(
            _models.import_model(spec, source, role, renorm_model=renorm_model))


# Old names (phase 6-quater, PR 6): see astrocook/_renamed.py
__getattr__ = renamed(__name__, {"RecipeFileV2": "RecipeFile"})
