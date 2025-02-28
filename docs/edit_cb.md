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

###  Import system list
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.import_systs</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>source</code>: Source session</li>
          <li><code>mode</code>: Mode (replace or append)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "import_systs",
  "params": {
    "source": "0",
    "mode": "'replace'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Import a system list into the currently selected session.*

This recipe imports a list of absorption systems from an open session `session` into the current one. Depending on `mode`, the list will either replace or be appended to an existing list of systems in the current session.


###  Import telluric model
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.import_telluric</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>source</code>: Source session</li>
          <li><code>col</code>: Telluric model column</li>
          <li><code>merge_cont</code>: Merge telluric model into continuum</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "import_telluric",
  "params": {
    "source": "0",
    "col": "'telluric_model'",
    "merge_cont": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Import a telluric model into the currently selected session and optionally merge it into the continuum.*

This recipe imports a telluric model from an open session `session` into the spectrum of the current session. The telluric model must be contained in a column `col` with the same length of the spectrum. If `merge_cont` is `True`, the model with be merged with the continuum of the spectrum and it will be used to normalize the flux when fitting absorption lines.


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
