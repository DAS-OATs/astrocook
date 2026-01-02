# Editing Data

Spectra rarely come ready-to-use out of the box. You might need to update metadata, rescale flux values, isolate specific regions, or combine multiple exposures. This tutorial covers the **Edit** menu, which serves as the "Swiss Army knife" for these operations.

(editing-properties)=
## 1. Setting Session Properties

Before running complex algorithms, ensure your session metadata is correct. Many Astrocook recipes (like the Ly-$\alpha$ forest analysis) rely on the emission redshift ($z_\mathrm{em}$) to determine rest-frame wavelengths.

1.  Select your session in the **Session List** (left sidebar).
2.  Go to **Edit > Set Properties...**.
3.  A dialog will appear. Here you can update:
    * **Name:** The label used in the list.
    * **Object:** The target name (often from the FITS header).
    * **Emission Redshift ($z_\mathrm{em}$):** Set this carefully!
    * **Resolution ($R$):** You can enter a resolving power (e.g., `50000`) or a pixel width (e.g., `3px`).

```{image} ../_static/editing_set_properties.png
:alt: The Set Properties dialog
:align: center
```

(editing-arithmetic)=
## 2. Arithmetic on Columns

Astrocook allows you to perform mathematical operations on your data columns (Flux `y`, Error `dy`, Wavelength `x`, etc.) using NumPy-style syntax.

### Modifying Flux

Suppose you want to rescale your flux by a factor of $10^{17}$ for easier readability:

1. Go to **Edit > Apply Expression...**.
2. **Target Column:** `y` (This is the column we will overwrite).
3. **Expression:** `y * 1e17` (This multiplies the current values by $10^{17}$).
4. Click **Run.**

### Creating New Columns

You can also create new auxiliary columns. For example, to create a column representing the Signal-to-Noise Ratio:

1. **Target Column**: `snr` (A new name).
2. **Expression:** `y / dy`.
3. Once created, you can visualize this new column using the Aux. Column dropdown in the [Plot Controls](#getting-started-customizing).

```{image} ../_static/editing_apply_expression.png
:alt: Applying a mathematical expression to a column
:align: center
```

## 3. Slicing Sessions (Split & Extract)

Sometimes you only want to work on a specific piece of a spectrum, such as the Ly-$\alpha$ forest or a specific metal line.

### Manual Split

To extract a custom range into a _new_ session:

1. Go to **Edit > Split Out Region...**.
2. Enter a boolean expression, such as `(x > 400) & (x < 500)` (assuming units of nm).
3. A new session containing only that data range will appear in your list.

### Quick Extraction (Presets)

Astrocook knows standard quasar regions. If your $z_\mathrm{em}$ is set:

1. Go to **Edit > Extract Preset Region...**.
2. Choose a preset like `lya_forest` (the region between Ly-$\beta$ and Ly-$\alpha$ emission).
3. The software automatically calculates the observed wavelength bounds based on $z_\mathrm{em}$ and creates a new session for you.

## 4. Combining Sessions (Stitch & Co-add)

If you have multiple exposures of the same object, you can combine them directly from the Session List.

### Co-adding Spectra

This operation defines common grid and rebins all spectra at once into it (weighting the contributions by their inverse variance and the amount of overlap with the final bins).

1. Load all the spectra you wish to combine.
2. In the **Session List**, hold `Ctrl` (or `Cmd` on Mac) and click to select multiple sessions.
3. Right-click on the selection.
4. Choose **Co-add N sessions...**.
5. A dialog will ask for the grid parameters (e.g., step size `dx`).
6. A new, combined session will be created.

### Stitching

If you have spectra covering _different_ wavelength ranges (e.g., Blue and Red arms of a spectrograph), use **Stitch**:

1. Select the sessions as above.
2. Right-click and choose **Stitch N sessions...**.
3. This concatenates the data arrays without resampling.

```{image} ../_static/editing_multiselect_coadd.png
:alt: Context menu for co-adding multiple sessions
:align: center
```