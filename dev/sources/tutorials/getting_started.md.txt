# Getting Started

Welcome to **Astrocook**! This guide will walk you through the graphical interface, from loading your first spectrum to performing basic analysis.

## 1. Launching the Application

Start Astrocook from your terminal:

```bash
astrocook
```

The Main Window will appear, divided into three sections:

1.  **Session List (Left):** Tracks the history of your modifications (Undo/Redo steps).
2.  **Plot Area (Center):** Interactive visualization of the spectrum.
3.  **Plot Controls (Right):** Toggles for errors, continuum, models, and redshift cursors.

## 2\. Loading Data

To load a spectrum, go to **File \> Open Spectrum...** or press `Ctrl+O`.
Astrocook supports:

  * **FITS files** (`.fits`)
  * **Astrocook Sessions** (`.acs`, `.acs2`)

:::{note}
When you load a file, Astrocook automatically creates a new **Session**. You can have multiple sessions open simultaneously (e.g., to compare different quasars).
:::

## 3\. The Plot Interface

Once data is loaded, the Right Sidebar becomes active. Here you can toggle visibility for different data components:

  * **1-sigma error:** Toggles the grey shading representing flux uncertainty.
  * **Continuum:** Toggles the dashed line for the fitted continuum.
  * **Absorption model:** Toggles the red line showing Voigt profile fits.
  * **Systems:** Toggles vertical markers for identified absorption systems.

### Navigation

  * **Pan:** Click and drag with the **Left Mouse Button**.
  * **Zoom:** Right-click and drag to define a zoom rectangle.
  * **Reset:** Click the **Home** icon in the toolbar to reset the view.

## 4\. Basic Analysis: Smoothing

Let's improve the signal-to-noise ratio visually by smoothing the spectrum.

1.  Go to **Flux \> Smooth Spectrum...**.
2.  A dialog will appear asking for parameters.
3.  **Sigma (km/s):** Enter `100.0`.
4.  Click **Run**.

Astrocook processes this "Recipe" in the background. Notice that a new entry appears in the **Session List** on the left.

:::{tip}
**Undo/Redo:**
If you don't like the result, simply press `Ctrl+Z` (Undo) or click the previous state in the Session List. Astrocook saves the full state of the spectrum at every step\!
:::

## 5\. Estimating the Continuum

To analyze absorption lines, we need to normalize the spectrum.

1.  Go to **Continuum \> Auto-estimate Continuum...**.
2.  Leave the default parameters (e.g., Kappa=2.0) and click **Run**.
3.  The blue line (Flux) will now be overlaid with a black dashed line (Continuum).

To view the normalized flux:

1.  Look at the **Plot Controls** (Right Sidebar).
2.  Check the box **Normalize F**.
3.  The Y-axis changes to `F / Continuum`.

## Advanced: Scripting API

Every action you take in the GUI corresponds to a command in the Astrocook API. You can automate your workflow using Python scripts.

:::{dropdown} Click to see the Python Code equivalent
The actions above can be replicated in a script using the `SessionV2` API:

```python
from astrocook import SessionV2

# 1. Load the session
sess = SessionV2.open_new(file_path="spec_qso.fits", name="QSO_1", gui_context=None, format_name="generic")

# 2. Smooth the flux (Flux Menu)
# This returns a NEW session object (the API is immutable)
sess = sess.flux.smooth(sigma_kms=100.0)

# 3. Auto-estimate Continuum (Continuum Menu)
sess = sess.continuum.estimate_auto(
    smooth_len_lya=5000.0, 
    kappa=2.0
)

# 4. Save the result
sess.save("output_session.acs2")