# Continuum Fitting

A robust continuum fit is the backbone of any absorption line analysis. Astrocook provides both automated algorithms for quick estimation and interactive tools for fine-tuning. This tutorial covers the **Continuum** menu and the manual editing tools.

## 1. Automatic Estimation

For most quasar spectra, the automated pipeline provides an excellent starting point. It uses an iterative kappa-sigma clipping algorithm to distinguish absorption lines from the continuum.

:::{note}
This recipe requires **Emission Redshift ($z_{em}$)** to be defined. If missing, Astrocook will automatically prompt you to set them via the [Set Properties](#editing-properties) dialog before proceeding.
:::

1. Go to **Continuum > Auto-estimate Continuum...**.
2. Adjust the parameters if necessary:
   - Smoothing (Ly-a / Out): Defines the scale of the smoothing spline. The default `5000` km/s for the Ly-alpha forest and `400` km/s for the red side usually work well.
   - Kappa: The sigma threshold for clipping absorption (default `2.0`). Lower values make the algorithm more aggressive at finding lines.
   - Fudge: Leave as `auto`. This calculates a correction factor to ensure the residuals of the unabsorbed regions are centered on zero.
3. Click **Run**.

Astrocook will generate a new continuum curve (black dashed line). If you have a `model` (absorption lines) already defined, it will automatically re-normalize it to match the new continuum.

:::{tip}
You can use the **Data Inspector** (Right-click on the plot) to examine the exact continuum values (`cont`) and residuals at any pixel. See [The Data Inspector](getting_started.md#data-inspector).
:::

```{image} ../_static/continuum_auto.png
:alt: Applying a mathematical expression to a column
:align: center
```

### Step-by-Step Method

If you need more control, you can run the process in two steps:

1. **Continuum > Find Absorbed Regions...**: Creates a mask (`abs_mask`) separating line/continuum.
2. **Continuum > Fit Continuum to Mask...**: Interpolates the continuum only over the unmasked regions. This allows you to manually edit the `abs_mask` (via [Apply Expression](#editing-arithmetic)) between steps if the automatic masking missed something.

(continuum-manual)=
## 2. Interactive Manual Editing

Automatic fits sometimes fail near complex emission lines or spectral edges. You can correct these manually using the **Right Sidebar**.

### Starting the Editor

1. Look at the **Edit Continuum** section in the right sidebar.
2. Click the **Start** button.
3. **Knots** (black dots) will appear along the current continuum.

:::{note}
If no continuum exists yet, Astrocook will prompt you to run the **Auto-estimate** routine first. You cannot edit the raw flux directly.
:::

:::{tip}
Manual editing modifies the session in place. If you want to compare your manual fit against the automatic one, use the **Duplicate** command in the Session List context menu to create a backup copy before you start editing.
:::

### Modifying Knots

- **Move**: Left-click and drag a knot to adjust the continuum level at that wavelength.

- **Add**: Right-click anywhere on the plot to add a new knot.

- **Remove**: Right-click on an existing knot to delete it.

### Controlling Stiffness (Slider & Spacing)

You can adjust the **knot density** using the slider or the **Spacing** text box.

- **Slider**: Drag left to increase density (closer knots) or right to decrease it (smoother curve).
- **Spacing Box**: Type a precise velocity spacing (in km/s) and press Enter.

:::{important}
Changing the spacing may affect the overall shape of the continuum. If you want to discard your manual moves and revert to the original continuum shape, click the **Reset** button (next to *Save*). This restores the knots to their default distribution based on the saved continuum.
:::

```{image} ../_static/continuum_manual_edit.png
:alt: Applying a mathematical expression to a column
:align: center
```

### Saving

When you are satisfied with the shape:

1. Click the **Save** button (which replaced "Start").
2. The new continuum is saved to the session, and the knots disappear.

## 3. Power-Law Fitting

For quasars, it is often useful to fit a power-law continuum ($F_\lambda\propto\lambda^{-\alpha}$) to estimate the underlying emission before determining the detailed shape.

:::{note}
This operation assumes your spectrum is **flux calibrated** (the overall shape is physically meaningful). It also requires **Emission Redshift ($z_{em}$)** to be defined. If missing, Astrocook will automatically prompt you to set them via the [Set Properties](#editing-properties) dialog before proceeding.
:::

1. Go to **Continuum > Fit Power-Law...**.
2. **Regions:** Enter a list of comma-separated wavelength ranges (in **rest-frame** nm) known to be free of emission lines.
   - Default: `128.0-129.0, 131.5-132.5, 134.5-136.0` (standard windows outside the Lyman forest).
3. Click **Run**.

This creates a new column `cont_pl` which is **automatically displayed** on the plot (via the *Aux. Column* selector). You can use it as a baseline for further fitting.