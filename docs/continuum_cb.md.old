---
layout: default
title: Continuum cookbook
parent: Cookbooks
nav_order: 2
math: mathjax2
---

# Continuum cookbook
{: .no_toc}

This cookbook contains utilities to detect absorption features in the spectrum and remove them, to determine the continuum level before absorption.


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Clip flux in spectrum
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.flux_clip</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>hwindow</code>: Half-window size in pixels for running mean</li>
          <li><code>kappa</code>: Number of standar deviations for clipping</li>
          <li><code>iter</code>: Number of iterations</li>
          <li><code>std</code>: Standard deviation for gaussian convolution (km/s)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "flux_clip",
  "params": {
    "hwindow": "100",
    "kappa": "5",
    "iter": "100",
    "std": "500"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Discriminate absorbed spectrum bins by applying a kappa-sigma clipping within a running window.*

The recipes computes the mean of `y` in a runniung window across the spectrum and saves it in column `y_rm`. It then computes the deviation $$\Delta y =$$ `y` $$-$$ `y_rm` for all spectral bins and finds the outliers in this distribution.

Since the recipe is meant to operate on quasar spectra, it assumes that the absorption features are narrow and widespread, while the emission feature are wide enough to be considered a part of the continuum. For this reason, the recipe is biased towards detecting outliers with a *negative* deviation (lower outliers, i.e. absorbed bins) instead of a *positive* deviation (upper outliers, most likely residuals of sky line or cosmic ray hit subtraction).

In practice, a bin is regarded as an outlier if it has either $$\Delta y < -$$`k` $$\times$$ `dy` or $$\Delta y > 2$$`k` $$\times$$ `dy`. Notice the two in the second formula, meaning that a bin must deviate twice as much in the upper direction than in the lower direction to be considered as an outlier.

The clipping is iterated several times, always removing the outliers and re-computing `y_rm`, until no more outliers are found or a maximum number of iterations `iter` is reached.

A column `y_abs` is created and set equal to `1` for all the lower outliers (absorbed bins), `0` elsewhere. Similarly, a column `y_em` is created for the upper outliers (“emitted” bins) and a column `y_cont` is created for the remaining bins (continuum bins). The latter is smoothed in the velocity space with a gaussian filter with standard deviation `std`, and the result is saved in column `cont` of the spectrum as the current estimate of the emission continuum.

###  Find lines
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.lines_find</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>std_start</code>: Start standard deviation of the gaussian (km/s)</li>
          <li><code>std_end</code>: End standard deviation of the gaussian (km/s)</li>
          <li><code>col</code>: Column to convolve</li>
          <li><code>kind</code>: Kind of extrema ('min' or 'max')</li>
          <li><code>kappa_peaks</code>: Number of standard deviations</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>append</code>: Append lines to existing line list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "lines_find",
  "params": {
    "std_start": "100.0",
    "std_end": "0.0",
    "col": "'y'",
    "kind": "'min'",
    "kappa_peaks": "5.0",
    "resol": "null",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Create a line list by convolving a spectrum with different gaussian profiles and finding the peaks in the convolution.*

The recipe is an extension of [peaks_find](#find-peaks), using convolution to detect peaks of different width. The gaussian profiles are applied in velocity space and used to filter out all fluctuations below a given scale.

The recipe proves most effective when a wide range of line widths is spanned, from the widest to the narrowest (typically, with gaussian standard deviation ranging from 100 to 0 km/s). A growing list of detected peaks is created and updated while spanning the range of line widths, eliminating duplicated findings.

The final list of `x`, `xmin`, `xmax`, `y`, and `dy` values for the peaks are saved in the current line list. If `append` is `True`, the new peaks are appended to the existing line list (if present), otherwise they replace it.

N.B. If the same line is found more than once, using different line widths for the gaussian filtering, it is saved just once in the line lists, and the value of `xmin` and `xmax` are the positions of the closest extrema (the closes maxima if the peak is a minimum, and vice-versa) when the line was first detected. This is the reason for spanning line widths in decreasing order: if a wide line is detected at first with a small value of gaussian standard deviation, it may be assigned a `xmin`-`xmax` range much smaller than the region it actually covers (because the filter is too fine to smooth out the noise around the peak), and the range is preserved even if the line is detected again with wider gaussians.

###  Continuum from nodes
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.nodes_cont</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>delta_x</code>: Size of slices (km/s)</li>
          <li><code>kappa_nodes</code>: Number of standard deviation away from the window average</li>
          <li><code>smooth</code>: Smoothing of the spline</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "nodes_cont",
  "params": {
    "delta_x": "500",
    "kappa_nodes": "5.0",
    "smooth": "0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Estimate a continuum by extracting, cleaning, and interpolating nodes from regions not affected by lines.*

This recipe combines [nodes_extract](#extract-nodes), [nodes_clean](#clean-nodes), and [nodes_interp](#interp-nodes) to estimate the continuum level once the lines have been detected.

###  Update line list
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.lines_update</code></td>
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
  "recipe": "lines_update",
  "params": {
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Update the line list with the parameters obtained by fitting the absorption systems.*

This recipe is meant to be run after the absorption systems have been fitted with Voigt profiles. It does not work if the list of systems comes from a different line lists.

For each systems, the best-fitting Voigt parameters are propagated to the line list and the FWHM of the lines is computed.

###  Find peaks
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.peaks_find</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>col</code>: Column where to look for peaks</li>
          <li><code>kind</code>: Kind of extrema ('min' or 'max')</li>
          <li><code>kappa</code>: Number of standard deviations</li>
          <li><code>append</code>: Append peaks to existing line list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "peaks_find",
  "params": {
    "col": "'conv'",
    "kind": "'min'",
    "kappa": "5.0",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Find the peaks in a spectrum column.*

Peaks are defined as the extrema (minima or maxima) of `y` column that are more prominent than `kappa` times the value of `dy` with respect to neighbouring bins.

The recipe finds all the peaks in the spectrum and saves their `x`, `xmin`, `xmax`, `y`, and `dy` values in the current line list. If the peak is a minimum, `xmin` and `xmax` are determined as the position of the closest maxima, and vice-versa if it is a maximum. If `append` is `True`, the new peaks are appended to the existing line list (if present), otherwise they replace it.

A column `line_mask` is created in the spectrum (if not present) and set to `1` within the `xmin`-`xmax` range of each detected peak.

###  Extract nodes
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.nodes_extract</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>delta_x</code>: Size of slices</li>
          <li><code>xunit</code>: Unit of wavelength or velocity</li>
          <li><code>mode</code>: Mode ('std' for creating nodes from spectrum, 'cont' for placing nodes across existing continuum)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "nodes_extract",
  "params": {
    "delta_x": "500",
    "xunit": "km / s",
    "mode": "'std'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Extract nodes from a spectrum, optionally masking lines.*

Depending on `mode`, the recipe can be used either to create a new set of nodes to estimate the continuum (`std`) or to place nodes across an existing continuum profile (`cont`). In the latter case, the spectrum must have a `cont` column beforehand.

In `std` mode, the recipe splits the spectrum in a number of slices and extracts a node for each slice by averaging the values of `x`, `y`, and `dy` within it. The node list is formatted like a normal spectrum.

Slices are defined in either wavelength or velocity space, as specified by the chosen `xunit`. They are equally spaced, with their size determined by `delta_x`.

If the spectrum has a `lines_mask` column, only the regions where this mask is `0` are considered when etracting the position of the nodes. This means that a slice may not be assigned a node if it is completely masked.

If the spectrum has a `deabs` column, containing the de-absorbed flux computed by fitting and removing the absorption features, this column is used instead of `y` to determine the `y` position of the nodes.

In `cont` mode, the `cont` column is used instead of `y` to determine the `y` position of the nodes, and the `lines_mask` is ignored.

###  Clean nodes
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.nodes_clean</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>kappa</code>: Number of standard deviation away from the window average</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "nodes_clean",
  "params": {
    "kappa": "5.0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Remove outliers from a list of nodes.*

The recipe applies the peak finding algorithm of [`peaks_find`](#find-peaks) to the list of nodes, to detect outliers and remove them iterately. Outliers are defined as the peaks that are more prominent than `kappa` times the value of `dy` with respect to the neighbouring nodes.

###  Interpolate nodes
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.nodes_interp</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>smooth</code>: Smoothing of the spline</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "nodes_interp",
  "params": {
    "smooth": "0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Interpolate nodes to estimate the emission level.*

The nodes in the current nodes list are interpolated with [`scipy.interpolate.UnivariateSpline`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html). The interpolating function is sampled at all values of `x` of the spectrum and the line list (if present), and is saved as the current continuum level in column `cont`.

###  Correct spectrum for Lyman-alpha opacity
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.lya_corr</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>zem</code>: Emisson redshift</li>
          <li><code>input_col</code>: Column to correct</li>
          <li><code>mode</code>: Correction mode ('basic' or 'inoue')</li>
          <li><code>logN_col</code>: Threshold for logarithmic column density</li>
          <li><code>percentile</code>: Percentile to compute the threshold from the column density distribution (only if logN_col is None)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "lya_corr",
  "params": {
    "zem": "zem",
    "input_col": "'y'",
    "mode": "'basic'",
    "logN_thres": "100",
    "percentile": "100"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Correct the spectrum flux for the effective Lyman-alpha opacity.*

Depending on `mode`, the correction factor is computed either from a simple integration of the H <span style="font-variant:small-caps;">i</span> growth function (`basic`), or using the prescriptions by Inoue et al. (2014) (`inoue`).

The `basic` mode require an upper integration limit that can be either provided through `logN_thres` or extracted from the distribution of H <span style="font-variant:small-caps;">i</span> column densities, from the table of fitted absorption systems (if present). In the latter case, the column density corresponding to a given `percentile` of the distribution is used as the upper integration limit.

The correction factor is multiplied to column `input_col` of the spectrum and the line list (if present), and the result is saved in column `[input_col]_taucorr`.

