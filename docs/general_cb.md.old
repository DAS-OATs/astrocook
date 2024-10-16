---
layout: default
title: General cookbook
nav_order: 0
math: mathjax2
---

# General cookbook
{: .no_toc}

This cookbook contains utilities to manipulate sessions, mask the spectrum, estimate spectral quality parameters, and perform basic operations like rebinning and convolution.


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Equalize two sessions
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.equalize</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>xmin</code>: Minimum wavelength (nm)</li>
          <li><code>xmax</code>: Maximum wavelength (nm)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "equalize",
  "params": {
    "xmin": "xmin",
    "xmax": "xmax",
    "_sel": "''"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Equalize the spectrum of two sessions, based on their flux ratio within a wavelength window.*

By default, the last-selected spectrum is equalized to the first-selected one (which is left unchanged). Equalization is done in place, without creating a new session.

To compute the rescaling factor, the recipe takes the medians of the `y` columns of the two spectra between `xmin` and `xmax`. The `y` and `dy` columns of the second spectrum are then multiplied by $$ \textrm{med}($$`y`$$_1)/\textrm{med}($$`y`$$_2)$$.

N.B. To select sessions, either click on the session window or provide a list through the hidden parameter `_sel`.

###  Combine two or more sessions
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.combine</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>name</code>: Name of the output session</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "combine",
  "params": {
    "name": "'*_combined'",
    "_sel": "''"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Create a new session combining the spectra from two or more other sessions.*

The recipe collects all the bins from the original spectra and puts them all together in the new spectrum. The bins retain their original size (defined by `xmin` and `xmax`), so they may overlap in the final spectrum. By default, they are ordered by ascending `x`.

All other structures from the original sessions (line lists, etc.) are not propagated to the new one.

N.B. To select sessions, either click on the session window or provide a list through the hidden parameter `_sel`.

###  Mask the spectrum
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.mask</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>col</code>: Column with the mask</li>
          <li><code>cond</code>: Condition</li>
          <li><code>new_sess</code>: Create a new session from masked spectrum</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "mask",
  "params": {
    "col": "'mask'",
    "cond": "''",
    "new_sess": "true",
    "masked_col": "'x'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Create a mask applying a specified condition to the spectrum bins.*

The expression in `cond` is translated into a boolean condition by the [`ast`](https://docs.python.org/3/library/ast.html) module. Expressions like `c>10` or `1500<c<2000` are supported, where `c` is a column of the spectrum.

The condition is checked on all spectrum bins and a new column `col` is populated with the results. No information is deleted from the input spectrum. If `new_sess` is `True`, a new session is created, containing a masked version of the input spectrum. In this masked spectrum, the column `y`, `dy`, and optionally `cont` are set to `numpy.nan` in all bins where the condition is false. All other structures from the original session (line lists, etc.) are not propagated to the new one.

###  Mask telluric absorption
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.telluric_mask</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>shift</code>: Shift to the heliocentric frame (km/s)</li>
          <li><code>apply</code>: Apply mask to flux</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "telluric_mask",
  "params": {
    "shift": "0",
    "apply": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Mask spectral regions affected by telluric absorptions.*

The regions were determined by Tobias M. Schmidt from ESPRESSO data and are saved in `telluric.dat`. They are resampled into the current `x` grid and used to populate a `telluric` column, which is set to `1` inside the regions and to `0` elsewhere.

If `apply` is `True`, `y` is set to `numpy.nan` in all bins where `telluric` is 1.

###  Estimate the SNR
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.snr_est</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        –
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "snr_est",
  "params": {
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Estimate the signal-to-noise ratio per pixel.*

A `snr` column is populated with `y`/`dy` ratios computed for all spectrum bins.

###  Estimate resolution
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.resol_est</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>px</code>: Number of bins per resolution element</li>
          <li><code>update</code>: Update column 'resol'</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "resol_est",
  "params": {
    "px": "3",
    "update": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Assign a resolution to spectral bins, assuming that the spectrum is designed to have a fixed number of bins per resolution element.*

This recipe is useful to populate the `resol` column in a spectrum (needed to fit the absorption systems) when it is empty, and information about the original sampling of the data is available. It does *not* try to infer the resolution from, e.g., the width of unresolved spectral feature.

###  Estimate error from RMS
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.rms_est</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>hwindow</code>: Half-size in pixels of the running window</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "rms_est",
  "params": {
    "hwindow": "100"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Estimate flux error by computing the root-mean-square (RMS) of the flux within a running window.*

The RMS is computed over `y` values and is saved in `y_rms`. It may be useful to compare the latter with `dy` to check that the formal error is consistent with the actual dispersion of `y` values.

###  Re-bin spectrum
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.rebin</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>xstart</code>: Start wavelength (nm)</li>
          <li><code>xend</code>: End wavelength (nm)</li>
          <li><code>dx</code>: Step in x</li>
          <li><code>xunit</code>: Unit of wavelength or velocity</li>
          <li><code>norm</code>: Return normalized spectrum, if continuum exists</li>
          <li><code>filling</code>: Value to fill region without data</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "rebin",
  "params": {
    "xstart": "null",
    "xend": "null",
    "dx": "10.0",
    "xunit": "km / s",
    "norm": "true",
    "filling": "nan"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Apply a new binning to a spectrum, with a constant bin size.*

The algorithm for re-binning is described in [Cupani et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016SPIE.9913E..1TC/abstract). It properly weights the flux contributions from the old bins to the new ones, also when the former overlap with each other (as it happens when several exposures of the same object are combined into a single spectrum).

The new grid is designed to fully cover the original range of the spectrum (when `xstart` and `xend` are `None`) or a specified range (useful when different spectra must be re-binned to the same grid). It is defined in either wavelength or velocity space, as specified by the chosen `xunit`. Any gap in the original binning are filled with a specified `filling` value, to ensure that the final grid is equally spaced.

Columns `y` and `dy` of the input spectrum are both re-binned to the new grid. If a column `cont` is present and `norm` is `True`, `y` and `dy` are normalized to `cont` in the re-binned spectrum.

A new session is created with the re-binned spectrum. All other structures from the original session (line lists, etc.) are not propagated to the new one.

###  Convolve with gaussian
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.gauss_convolve</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>std</code>: Standard deviation of the gaussian (km/s)</li>
          <li><code>input_col</code>: Input column</li>
          <li><code>output_col</code>: Output column</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "gauss_convolve",
  "params": {
    "std": "20.0",
    "input_col": "'y'",
    "output_col": "'conv'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Convolve a spectrum column with a gaussian profile.*

The convolution is computed in velocity space, using the Fast Fourier Transform.

###  Compute the CCF
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.flux_ccf</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>col1</code>: First column</li>
          <li><code>col2</code>: Second column</li>
          <li><code>dcol1</code>: Error for first column</li>
          <li><code>dcol2</code>: Error for second column</li>
          <li><code>vstart</code>: Start velocity</li>
          <li><code>vend</code>: End velocity</li>
          <li><code>dv</code>: Velocity step</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "flux_ccf",
  "params": {
    "col1": "'y'",
    "col2": "'y'",
    "dcol1": "'dy'",
    "dcol2": "'dy'",
    "vstart": "-20",
    "vend": "20",
    "dv": "0.1"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Convolve the cross-correlation function (CCF) between two spectrum columns.*

The recipe is designed to work on flux densities. Typically, the first column is `y` and the second column contains the flux density from a different spectrum with the same wavelength binning. The second columns can also be `y`: in this case, the recipe computes the auto-correlation instead of the cross-correlation.

The CCF is computed in velocity space, shifting `col2` with respect to `col1` within the range `vstart`-`vend` and with step @dv. The columns are resampled while shifting, to accomodate for values of @dv much smaller than the spectrum bin size.

The CCF is saved in a NumPy binary file `SESS_ccf.npy`, with `SESS` the name of the session.

###  Compute statistics of the CCF
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.flux_ccf_stats</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>n</code>: Number of realizations</li>
          <li><code>col1</code>: First column</li>
          <li><code>col2</code>: Second column</li>
          <li><code>dcol1</code>: Error for first column</li>
          <li><code>dcol2</code>: Error for second column</li>
          <li><code>vstart</code>: Start velocity (km/s)</li>
          <li><code>vend</code>: End velocity (km/s)</li>
          <li><code>dv</code>: Velocity step (km/s)</li>
          <li><code>fit_hw</code>: Half-window used for fitting the CCF (km/s)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "flux_ccf_stats",
  "params": {
    "n": "10.0",
    "col1": "'y'",
    "col2": "'y'",
    "dcol1": "'dy'",
    "dcol2": "'dy'",
    "vstart": "-20",
    "vend": "20",
    "dv": "0.1",
    "fit_hw": "1.0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Compute statistics for the peak of the cross-correlation function (CCF) by bootstrapping a number of realizations for the spectrum.*

Realizations are created by selecting entries at random, preserving wavelength order and rejecting duplicates (compare with Peterson et al. 1998).

The recipe computes the CCF between the original flux and the flux of each realization. A gaussian is fit to the CCF within a window around 0 (in velocity space) to determine the position of the peak. The distribution of peak positions is saved in a NumPy binary file `SESS_ccf_stats.npy`, with `SESS` the name of the session.
