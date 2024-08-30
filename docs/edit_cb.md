---
layout: default
title: Edit cookbook
parent: Cookbooks
nav_order: 0
---

# Edit cookbook
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Modify columns

🚧

###  Import system lists

🚧

### Extract

🚧

###  Mask
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.mask</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>shift</code>: Shift to the barycentric frame (km/s)</li>
          <li><code>tell</code>: Mask telluric lines</li>
          <li><code>sky</code>: Mask sky lines</li>
          <li><code>cond</code>: Condition</li>
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
    "shift": "0",
    "tell": "true",
    "sky": "true",
    "cond": "''"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Mask telluric lines, sky lines, or spectral regions defined by a condition.*

The user can mask telluric absorption lines and sky emission lines using `tell` and `sky` respectively, and using `shift` to specify by how much these lines have to be shifted with respect to the barycentric reference frame in which the spectrum is shown (this is the value of the barycentric correction frequently indicated in the spectrum header, expressed in km/s). The user can also mask other spectral regions with `condition`: in this case, the condition must be expressed in terms of one of the [spectrum table](tables.md#spectrum-table) columns (e.g. `x>500` to mask all wavelengths greater than 500 nm, `y<0.1` to mask all flux values less than 0.1 in the current flux units, etc.).

When a mask is applied, a new column is created (`telluric`, `sky`, or `mask`, respectively) and filled with 1 outside the masked regions and with 0 inside. Masking is destructive: where the mask column is 0, the flux values (column `y`) are changed into NaNs. This behaviour may be improved in a later version of the code. 

### Toggle log axes

🚧
