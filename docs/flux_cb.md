---
layout: default
title: Flux cookbook
parent: Cookbooks
nav_order: 1
math: mathjax2
---

# Flux cookbook
{: .no_toc}

Cookbook of utilities for flux calibration


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Scale y axis
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookFlux.y_scale</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>fact</code>: Multiplicative factor</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "y_scale",
  "params": {
    "fact": "1.0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Scale the y axis by a constant factor.*

The `y` and `dy` columns of the spectrum and the line list (if present) are multiplied by `fact`.

The scaling is done in place, without creating a new session.

###  Scale y axis by median
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookFlux.y_scale_med</code></td>
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
  "recipe": "y_scale_med",
  "params": {
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Scale the y axis by its median.*

The `y` and `dy` columns of the spectrum and the line list (if present) are multiplied by the median of the spectrum `y`.

The scaling is done in place, without creating a new session.

###  Scale y axis by its value at a given point
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookFlux.y_scale_x</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>x</code>: x (nm)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "y_scale_x",
  "params": {
    "x": "x"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Scale the y axis by its value at a given point.*

The `y` and `dy` columns of the spectrum and the line list (if present) are multiplied by the value of the spectrum `y` at a given `x`, computed with [`numpy.interp`](https://numpy.org/doc/stable/reference/generated/numpy.interp.html?highlight=interp#numpy.interp).

The scaling is done in place, without creating a new session.

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



