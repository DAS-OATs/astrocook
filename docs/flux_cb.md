---
layout: default
title: Flux cookbook
parent: Cookbooks
nav_order: 1
math: mathjax2
---

# Flux cookbook
{: .no_toc}

Cookbook of utilities to manipulate flux


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

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

🚧

###  Smooth spectrum

🚧

###  Rescale spectrum

🚧

###  De-redden spectrum
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookFlux.deredden</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>ebv</code>: Color excess</li>
          <li><code>rv</code>: Ratio of total selective extinction</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "deredden",
  "params": {
    "ebv": "0.03",
    "rv": "3.1"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Correct the spectrum flux for reddening due to extinction.*

The extinction is modeled with the parametrization of O'Donnell (1994), depending on the spectrum color excess $$E(B-V)$$ and ratio of total selective extinction $$R(V)=A(V)/E(B-V)$$. Column `y` of the spectrum is updated with de-reddened values.

🚧

###  Adjust magnitudes

🚧

###  Estimate SNR
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

###  Estimate RMS
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
