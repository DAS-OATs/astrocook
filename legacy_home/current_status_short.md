### Context: Astrocook V2 (Bug Fixing Mode)
* **Stack:** Python 3.10+, PySide6, Matplotlib, SciPy, Astropy.
* **Architecture:**
    * **Core:** Uses immutable data models (`SpectrumDataV2`, `SystemListDataV2`).
    * **Session:** State is managed by `SessionV2` (orchestrator) and `SessionHistory` (undo/redo).
    * **UI:** `MainWindowV2` holds a `central_stack` and floating sidebars.
* **Goal:**