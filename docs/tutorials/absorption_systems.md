# Absorption Systems

Detecting and fitting absorption lines is the core mission of Astrocook. The software provides a tiered approach: fully automated pipelines for initial detection, and a dedicated **System Inspector** for detailed Voigt profile analysis and manual refinement.

## 1. Automated Detection Pipelines

For most science cases, you should start with one of the automated pipelines in the **Absorbers** menu. These recipes scan the spectrum, identify features, and perform an initial fit.

### Ly-alpha Forest
If you are analyzing a high-redshift quasar ($z > 2$):
1.  Ensure your [Emission Redshift is set](#editing-properties) correctly.
2.  Go to **Absorbers > Auto-fit Ly-alpha Forest...**.
3.  **Threshold:** The default score (`0.1`) is usually fine.
4.  **Min b:** The minimum Doppler parameter (e.g., `10` km/s) helps filter out noise spikes or unphysically narrow lines (interlopers).
5.  Click **Run**.

This pipeline focuses on the region between the Ly$\beta$ and Ly$\alpha$ emission lines.

### Metal Doublets
For finding distinct metal systems (CIV, SiIV, MgII) redward of the Ly$\alpha$ emission:
1.  Go to **Absorbers > Auto-detect Doublets...**.
2.  **Multiplets:** You can customize the list (e.g., `CIV, SiIV, MgII`).
3.  Click **Run**.

This recipe searches for pairs of lines with the correct velocity separation and optical depth ratios.

### General Detection
For a more generic search (or if you are looking for specific ions):
1.  Go to **Absorbers > Identify Absorption Lines...**.
2.  This generates a list of candidate regions but *does not* fit them immediately.
3.  You can then visualize these regions by selecting `abs_ids` in the **Aux. Column** dropdown of the Plot Controls.

```{image} ../_static/absorbers_auto_pipeline.png
:alt: The Auto-fit Ly-alpha Forest dialog
:align: center
```

## 2. The System Inspector

Once you have systems (either from a pipeline or added manually), the **System Inspector** is your cockpit for detailed analysis.

To open it: **View > View System Inspector**.

The window is divided into two parts:

- **Left (System List)**: A table of all absorption components.
- **Right (Velocity Plot)**: A zoomed-in view of the selected system in velocity space.

### Navigating Systems

- **Select a row** in the table to jump to that component.
- The Velocity Plot will automatically center on the selected line (e.g., `CIV_1548`) and show all associated transitions (e.g., `CIV_1550`) stacked vertically.
- **Group View**: Check the "Group View" box above the table to filter the list. This hides unrelated systems and shows only the components physically linked to your selection (e.g., all lines in the same doublet).

```{image} ../_static/absorbers_system_inspector.png
:alt: The System Inspector showing a metal doublet
:align: center
```

## 3. Interactive Fitting & Editing

The System Inspector allows you to modify the model dynamically.

### Adding Components

If the auto-finder missed a line:

1. In the Velocity Plot (right panel), **Right-click** on the background where the line should be.
2. Select **Add [Ion] at z=...**.
3. A new component is added and fitted immediately.

:::{tip}
If you see multiple transitions in the menu (e.g., `Add CIV-SiIV`), this will add linked components for both ions at that redshift.
:::

### Modifying Parameters

You can edit values directly in the table:

- **Double-click** a cell (`z`, `logN`, `b`) to type a new value.
- The model (red line) updates instantly to reflect your change.

### Constraints (Fix/Free)

To fix a parameter during fitting:

1. Select the component(s) in the table.
2. **Right-click** to open the context menu.
3. Select **Freeze 'b'** (or `z`, `logN`).
4. The value will turn _italic_, indicating it is locked.

### Parameter Linking

Astrocook supports physical linking of parameters (e.g., tying the redshift of a CIV doublet).

1. Select **two** components in the table (hold `Ctrl` or `Cmd`).
2. Right-click on the parameter you want to link (e.g., the `z` column).
3. Select **Link z to [Other Component] (Value)**.
4. The value will turn **bold**, indicating it is dynamically calculated from the other component.

### Refitting

After manual edits, you often need to re-optimize the fit:

- **Refit Selected:** Right-click the table and choose "Refit Selected" to optimize only the highlighted lines.
- **Refit All:** Go to **Absorbers > Refit All Systems...** (in the main window) to re-optimize the entire spectrum, handling blends automatically.

```{image} ../_static/absorbers_linking.png
:alt: Linking parameters via the context menu
:align: center
```