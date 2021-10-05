---
layout: default
title: General cookbook
parent: Cookbooks
nav_order: 0
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

_Equalize the flux level of one session to another one_. You can select the two sessions clicking on the Sessions window or providing a list through the `xmax` is used to compute the median. Note that the first-selected session is left unchanged, while the other one is rescaled.

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

_Combine the spectra from two or more sessions_. You can select sessions clicking on the Sessions window or providing a list through the `_sel` parameter. A new session is created, with a new spectrum containing all entries from the spectra of the combined sessions. Other objects from the sessions (line lists, etc.) are discarded.

###  Create a spectral mask
        
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

_Create a spectral mask by applying a given condition_. The condition must be parsable by AST, with spectrum columns denoted by their names (e.g. 'x>400'). Optionally, a new session is created with the masked spectrum. Other objects from the old session (line lists, etc.) are discarded.

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

_Mask telluric absorption_

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

_Estimate the signal-to-noise ratio per pixel_.

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
          <li><code>px</code>: Number of pixels</li>
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

_Estimate spectral resolution assuming the spectrum has a fixed number of pixels per resolution element_.

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
          <li><code>hwindow</code>: Half-window size in pixels for running mean</li>
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

_Estimate flux error by computing the running RMS of the flux_.

###  Rebin spectrum
        
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

_Rebin a spectrum with a given step_. The step can be expressed in any unit of wavelength or velocity. Start and end wavelength may be specified, e.g. to align the rebinned spectrum to other spectra. If start or end wavelength are None, rebinning is performed from the first to the last wavelength of the input spectrum. A new session is created with the rebinned spectrum. Other objects from the old session (line lists, etc.) are discarded.

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

_Convolve a spectrum column with a gaussian profile using FFT transform_.

