---
layout: default
title: Continuum cookbook
parent: Cookbooks
nav_order: 2
---

# Continuum cookbook
{: .no_toc}

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

Select flux pixels to find outliers

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

Create a line list by convolving a spectrum with different gaussian profiles and finding the peaks in the convolved spectrum

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

Estimate a continuum by extracting, cleaning, and interpolating nodes from regions not affected by lines

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

Update line list after the systems have been fitted, copying fitting parameters and computing the line FWHM. The recipe only works if the systems were extracted from the line list that is to be updated.

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

Find the peaks in a spectrum column. Peaks are the extrema (minima or maxima) that are more prominent than a given number of standard deviations. They are saved as a list of lines.

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
          <li><code>mode</code>: Mode ('std' for extracting nodes from spectrum, 'cont' for converting continuum into nodes)</li>
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

Extract nodes from a spectrum. Nodes are averages of x and y in slices, computed after masking lines.

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

Clean the list of nodes from outliers.

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

Interpolate nodes with a univariate spline to estimate the emission level.

###  Correct flux for Lyman-alpha opacity
        
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

Correct flux for Lyman-alpha opacity, using the prescriptions by Inoue et al. 2014

###  Continuum from absorbers
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.abs_cont</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>zem</code>: Emisson redshift</li>
          <li><code>std</code>: Standard deviation of the gaussian (km/s)</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>mode</code>: Correction mode ('basic' or 'inoue')</li>
          <li><code>reest_n</code>: Number of re-estimation cycles</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "abs_cont",
  "params": {
    "zem": "zem",
    "std": "1000.0",
    "resol": "null",
    "mode": "'basic'",
    "reest_n": "4",
    "_refit_n": "0",
    "_percentile": "100",
    "_print_stats": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Estimate a continuum by iteratively fitting and removing absorbers

