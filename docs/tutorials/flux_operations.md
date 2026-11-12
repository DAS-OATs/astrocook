# Flux Operations

Once you've loaded your data and set the basic properties, you might need to process the flux arrays to increase signal-to-noise ratio (SNR) or prepare for multi-spectrum analysis. This tutorial covers the **Flux** menu operations.

## 1. Calculating Running Statistics

Sometimes you need to estimate the noise level directly from the data itself (e.g., if the error column is missing or unreliable).

1. Go to **Flux > Calculate Running StdDev...**.
2. **Input Column:** Usually `F` (the flux, the default).
3. **Window Size:** The number of pixels to use for the sliding window (e.g., `21`). The window is an odd number of pixels: an even one is increased by 1, and one below 3 is 3.
4. **Output Column:** Name for the new error column (default: `running_std`). An existing auxiliary column of that name is replaced; the columns of every spectrum (`λ`, `λmin`, `λmax`, `F`, `dF`) cannot be used.

This creates a new auxiliary column containing the local Root-Mean-Square (RMS) variations, which you can then visualize or use as a mask.

```{image} ../_static/flux_running_std.png
:alt: The running standard deviation of the flux
:align: center
```

## 2. Smoothing

Smoothing is the simplest way to visualize faint features by reducing high-frequency noise. Astrocook uses a Gaussian kernel for this operation.

1. Go to **Flux > Smooth Spectrum...**.
2. **Sigma (km/s)**: Enter the width of the Gaussian kernel.
   - A value of `100` km/s is often a good starting point for visualization.
   - Smaller values (e.g., `10-30` km/s) preserve narrow absorption lines.
3. Click **Run**.

The operation creates a _new_ session state with the smoothed flux. The original data is preserved in the history (undo with `Ctrl+Z` / `Cmd+Z`). Note that the error array (`dy`) is also smoothed to maintain consistency, and so are the continuum (`cont`), the model (`model`), the power-law continuum (`cont_pl`) and the telluric transmission (`telluric_model`), if present.

```{image} ../_static/flux_smooth.png
:alt: A smoothed spectrum
:align: center
```

## 3. Rebinning

Rebinning is more rigorous than smoothing. It combines adjacent pixels into larger bins: the flux of each new bin is the inverse-variance weighted mean of the old pixels it covers, each weighted by the fraction of it that falls in the bin. This is essential when you want to improve SNR for quantitative measurements or reduce the data size.

1. Go to **Flux > Rebin Spectrum...**.
2. **Step Size (dx)**: Enter the desired size of the new bins.
3. **Unit:** Choose the unit for `dx`: a velocity (`km/s`, a grid uniform in log wavelength) or a wavelength (`nm`).
4. **Wavelength Range (Optional)**: You can define a specific start (`xstart`) and end (`xend`) wavelength, in nm; the start must be below the end. Leaving them as None will use the full existing range.
5. **Kappa (Optional)**: Sigma-clipping threshold. Pixels more than `kappa` sigma *above* a running median of the input pixels (spikes, cosmic rays) are rejected before rebinning; absorption features are never clipped, and a bin is never left empty. The default is 5; `None` turns the clipping off.
6. Click **Run**.

The error of each new bin is the error of the weighted mean of the pixels it contains (it decreases as $\sigma/\sqrt{N}$ for $N$ equal pixels).

:::{note}
Unlike smoothing, rebinning changes the X-axis grid. The new flux values are calculated to preserve the physical energy/counts per bin.
:::

## 4. Equalizing

Sometimes the flux calibration between two exposures or arms of a spectrograph differs. You can scale the current session's flux to match a reference session.

1.  Go to **Flux > Equalize to Reference...**.
2.  **Reference Session:** Click the session you want to match in the list of the open sessions (or type its name); it is required.
3.  **Order:** Choose the complexity of the scaling:
    * `0` (Scalar): Multiplies flux by a constant factor.
    * `-1` (Spline): Fits a smooth curve to match the shape (the lowest value allowed).
    * `1` or higher: Fits a polynomial of that order.

:::{note}
Here `-1` means a spline; in **Co-add** (`equalize_order`) it means no scaling.
:::

## 5. Resampling

Resampling is different from rebinning. It interpolates your spectrum onto an arbitrary grid defined by _another_ open session. This is a critical prerequisite for comparing two spectra pixel-by-pixel (e.g., for creating a ratio or difference spectrum).

**Scenario:** You have a high-resolution spectrum (Session A) and a low-resolution standard (Session B), and you want to divide A by B.

1. Load both sessions.
2. Select Session A (the one you want to change).
3. Go to **Flux > Resample on Grid...**.
4. **Target Session:** Click "Session B" in the list of the open sessions (or type its name); it is required.
5. Click **Run**.

Session A will now have exactly the same wavelength points as Session B, allowing you to use [Apply Expression](#editing-arithmetic) to compute `F / SessionB.F`.

(flux-calibration)=
## 6. Flux Calibration

If your spectrum is not flux-calibrated (or the calibration is off), you can warp it to match known photometric data.

1.  Go to **Flux > Flux Calibrate...**.
2.  **Magnitudes:** Enter a comma-separated list of `Filter=Magnitude` pairs (AB magnitudes).
    * Example: `SDSS_g=19.2, SDSS_r=18.8`
    * The filters: `SDSS_u`, `SDSS_g`, `SDSS_r`, `SDSS_i`, `SDSS_z`, `Johnson_U`, `Johnson_B`, `Johnson_V`, `Johnson_R`, `Johnson_I`. A filter not in the list, or a pair not written as `Filter=Magnitude`, is refused before the recipe runs.
3.  Click **Run**.

Astrocook calculates a smooth correction curve to force the synthetic photometry of your spectrum to match the input magnitudes. With one magnitude the curve is a constant; with two or more it is linear between the central wavelengths of the filters and constant beyond them, and its values there are chosen so that every magnitude is matched exactly. If the magnitudes cannot be matched together by a positive curve (e.g. two strongly overlapping filters with very different magnitudes, or the same filter twice), a warning is logged and each filter's own factor is used, so the magnitudes are only approximately matched. The curve is applied to the flux, its error, the continuum, the model and the power-law continuum, which keep the unit of the flux (e.g. erg cm⁻² s⁻¹ nm⁻¹); the telluric transmission is dimensionless and is left unchanged. A flux in arbitrary units (see [the unit of the flux](#getting-started-flux-unit)) is read as if it were in erg cm⁻² s⁻¹ Å⁻¹, with a warning in the log, and takes that unit.

:::{tip}
After calibration, check the **Show iso-magnitude grid** box in the [Plot Controls](#getting-started-customizing) to verify the result against reference magnitude curves.
:::