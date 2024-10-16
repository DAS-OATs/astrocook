---
layout: default
title: Absorbers cookbook
parent: Cookbooks
nav_order: 3
---

# Absorbers cookbook
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Find lines
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.find_lines</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>kind</code>: Kind</li>
          <li><code>prominence</code>: Prominence</li>
          <li><code>append</code>: Append to existing line list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "find_lines",
  "params": {
    "kind": "'abs'",
    "prominence": "null",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Find absorption or emission lines, based on their prominence.*

Lines are found in column `y` of the spectrum. `kind` is either `abs` or `em`, to find absorption or emission lines respectively. `prominence` is defined as in `scipy.signal.find_peaks`; if the user does not specify it, it is defined as 5 times the value in column `dy` of the spectrum.

###  Model Ly-a forest

🚧

###  Model metals
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.model_metals</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>trans</code>: Transitions</li>
          <li><code>zem</code>: Emission redshift</li>
          <li><code>no_ly</code>: Exclude Lyman forest</li>
          <li><code>use_lines</code>: Use line list to define absorbers @url absorbers_cb.html#model-metals</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "model_metals",
  "params": {
    "series": "series",
    "zem": "zem",
    "no_ly": "true",
    "use_lines": "false"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

*Model metal absorbers, based on transition and emission redshift.*

The spectrum is scanned to identify absorption systems and fit them with Voigt profiles. If `series` is a doublet or multiplet (recommended), e.g. `CIV` or `MgII_2796,MgII_2803`, only absorbers with the right wavelength ratio(s) are modelled as bona fide systems and saved in the system table with their fitting parameters. Only wavelengths up to the emission redshift `zem` of the chosen `series` are considered. The user may decide to exclude the Lyman forest with `no_ly`. If `use_lines` is `True`, absorption systems are identified on an existing line list, otherwise they are identified by cross-correlating the flux spectrum with itself, after applying a shift corresponding to the wavelength ratio(s) under consideration.

When the system table is empty, the recipe tries to identify systems using the column `y` of the spectrum, otherwise it tries on the column `deabs` (to avoid considering absorbers already identified and modeled). The recommended practice is to use this recipe to identify the most common doublets, and to identify the remaining lines using [`identify_unknown`](absorbers_cb.md#identify-unknown-lines) instead. If the Lyman forest has already been modeled, the user may want to set `no_ly` to `False` and tentatively identify metals within the forest too.

###  Identify unknown lines

🚧

###  Check system list

🚧
