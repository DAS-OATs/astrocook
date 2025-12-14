### **Current Project Status: Astrocook V2 (Batch-Ready Beta)**

**Project Context:** Astrocook V2 is now a robust, scriptable analysis package (`astrocook/`) capable of both GUI interactivity and high-volume batch processing. V2 is the active development branch (`develop`), now verified against synthetic "stress tests."

**Current Focus:** **Scientific Validation & High-Volume Deployment** (Core architecture and fitting logic are stable).

### **I. Architecture & Directory Structure (Stable)**

* **Root Package:** `astrocook/`
* **Core (`astrocook/core/`):**
  * **Data:** Immutable models `SpectrumDataV2`, `SystemListDataV2`.
  * **Logic:** Orchestrator `SessionV2`, logic provider `SystemListV2`.
  * **State Management:** `SystemListV2` now includes a `fitting_context` manager to safely handle temporary fitting masks.
* **I/O Layer (`astrocook/io/`):**
  * **Registry Pattern:** `@register_loader` supports native formats.
  * **Legacy Bridge:** `adapter.py` supports legacy file migration.
  * **Robust Persistence:** The system now correctly serializes and deserializes fit metadata (`chi2`, `resol`) via the `.acs2` archive format using "comment smuggling" to bypass legacy constraints.
* **Scripting:** `batch_driver.py` (Controller) and `parallel_orchestrator.py` (Executor) enable headless analysis of 10k+ sightlines.

### **II. Physics & Fitting (`astrocook/fitting/`)**

* **`VoigtFitterV2`:** SciPy-based fitting engine.
  * **Smart Resolution:** Automatically detects if the `resol` column is in $R$ (Resolving Power) or km/s (FWHM) based on magnitude heuristics.
  * **SSOT Convolution:** A unified `convolve_flux` function ensures the Model, Plot, and Cursor all use identical convolution logic.
* **Optimization Logic (`SystemListV2.optimize_hierarchy`):**
  * **Fixed-Window AIC:** Calculates AIC on a static pixel mask to prevent "shifting goalpost" biases during iterative component addition.
  * **Taboo Search:** Uses a `rejected_mask` to prevent the optimizer from looping on failed candidates.
  * **Proximity Checks:** Prevents singular fits by rejecting candidates too close (`min_dv`) to existing lines.

### **III. GUI Implementation (`astrocook/gui/`)**

#### **1. Main Window (`main_window.py`)**
* **Stack:** PySide6 + Matplotlib.
* **Interaction:** Threaded recipe execution.

#### **2. System Inspector (`system_inspector.py`)**
* **Visuals:**
  * **Accurate Cursor:** The interactive Voigt cursor now correctly interpolates variable resolution columns to match the static model.
  * **Decluttered Plot:** Resolution labels removed; metadata moved to table tooltips.
* **Controls:**
  * **Table:** Tooltips now correctly display `Chi2` and `Resol` loaded from batch results (no longer "N/A").

### **IV. Recipes (`astrocook/recipes/`)**

Astrocook possesses and API layer exposing all recipes to the user, for scripting and LLM integration. Some examples:
* **`absorbers.py`:**
  * **`optimize_system`:** Wraps the core `optimize_hierarchy` logic, ensuring the visual model is updated (unmasked) after the fit completes.
* **`edit.py`:**
  * **`trim_common`:** New recipe to trim spectra to the velocity intersection of multiple chunks (HI + Metal 1 + Metal 2).

### **V. Recent Accomplishments**

1.  **Batch Architecture:** Successfully implemented `parallel_orchestrator` to process 10,000+ sightlines with robust logging suppression.
2.  **Scientific Robustness:** Verified `optimize_hierarchy` against synthetic "Train Wreck" data (SNR=10, 10 blended components); achieved ~80% recovery with no overfitting.
3.  **Persistence Fixes:** Solved the "N/A" metadata bug by updating `adapter.py` and `system_list_migration.py` to persist fit statistics in legacy tables.
4.  **UX Stability:** Fixed the vanishing cursor and synchronized resolution display between plots and tables.
5.  **Co-addition:** Solved critical aliasing ("jaggedness") in echelle co-addition by implementing a physics-aware pipeline that preserves 2D pixel geometry and uses vectorized drizzling.

### **VI. Task for this conversation**



### **VII. Files for Context**

