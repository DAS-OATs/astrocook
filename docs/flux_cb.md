---
layout: default
title: Flux cookbook
parent: Cookbooks
nav_order: 1
---

# Flux cookbook
{: .no_toc}

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

Scale the y axis by a constant factor. The spectrum and the line list are rescaled in place, without starting a new session.

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

Scale the y axis by its median. The spectrum and the line list are rescaled in place, without starting a new session.

###  Scale y axis by its value at a given x
        
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

Scale the y axis by its value at a given x.

###  Deredden spectrum
        
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
          <li><code>ebv</code>: Color excess E(B-V)</li>
          <li><code>rv</code>: Ratio of total selective extinction R(V)=A(V)/E(B-V)</li>
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

Deredden the spectrum using the parametrization by Cardelli, Clayton, and Mathis (1989) and O'Donnell (1994).

