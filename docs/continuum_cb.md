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

###  Clip flux
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookContinuum.clip_flux</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>ran</code>: Wavelength range (nm)</li>
          <li><code>smooth_len</code>: Smoothing length (km/s)</li>
          <li><code>kappa</code>: Number of sigma to reject absorber</li>
          <li><code>fudge</code>: Fudge factor to scale the continuum</li>
          <li><code>knots_dist</code>: Distance between knots (km/s)</li>
          <li><code>mode</code>: Update or replace</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "clip_flux",
  "params": {
    "ran": "'all'",
    "smooth_len": "400",
    "kappa": "2",
    "fudge": "'auto'",
    "knots_dist": "2000",
    "mode": "'update'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Estimate the continuum by clipping absorbers.*

Continuum is estimated by computing a running mean on `y` and iteratively rejecting negative outliers, defined by the following condition:

```math
y_\mathrm{running median} - y_\mathrm{outlier} > \kappa\times dy.
```

The continuum is further adjusted by a `fudge` factor which is either provided by the user or computed by minimizing the cumulative residuals in a running window. Finally, the continuum is smoothed with a gaussian filter of length `smooth_len` and saved in column `cont` of the spectrum.

Decrease `kappa` if the continuum follows too much the absorption systems. Decrease `smooth_len` if the continuum cuts out the emission lines. The continuum can be manually refined by adding/removing knots on the graph. `knots_dist` controls the distance of knots that are placed automatically over the continuum curve.

You can use `ran` to limit the continuum estimation to a specific wavelength range. If `mode` is `update`, the new continuum will be merged with the existing estimate (if present). If `mode` is `replace`, the existing estimate is deleted.

###  Fit power law

🚧

###  Correct for Ly-a opacity

🚧
