---
layout: default
title: General cookbook
parent: Cookbooks
nav_order: 1
---

# General cookbook
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

## Import structure

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>GUIPanelSession.struct_import</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>struct</code>: Structure</li>
          <li><code>mode</code>: Mode (<code>replace</code> or <code>append</code>)</li>
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Import a data structure from a session into the current one.** The structure to be imported is described by a string with the session number and the structure tag (`spec`, `lines`, `systs`) separated by a comma (e.g. `0,spec`, meaning "spectrum from session 0"). The imported structure is either replaced or appended to the corresponding one in the current session.

## Modify structures

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>GUIPanelSession.struct_modify</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>col_A</code>: Column A</li>
          <li><code>col_B</code>: Column B or scalar</li>
          <li><code>struct_out</code>: Output column</li>
          <li><code>op</code>: Binary operator</li>
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Modify a data structure using a binary operator.** An output column is computed applying a binary operator to two input columns, or an input column and a scalar. Columns are described by a string with the session number, the structure tag (`spec`, `lines`, `systs`), and the column name separated by a comma (e.g. `0,spec,x`, meaning "column x of spectrum from session 0").  They can be from different data structures only if they have the same length. If the output column already exists, it is overwritten.

## Extract region

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.region_extract</code></td>
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
  </tbody>
</table>

**Extract a spectral region.** The region between a minimum and a maximum wavelength is extracted from the data structures in the current session (these include the selected spectral range with all the lines and the absorption systems that fall within). A new session with the extracted data structures is created.

## Convert x axis

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.x_convert</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>zem</code>: Emission redshift</li>
          <li><code>xunit</code>: Unit of wavelength or velocity</li>
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Convert the x axis to wavelength or velocity units.** The x axis can be converted to any unit of wavelength or velocity (default: nm and km/s). When converting to and from velocity units the zero point is set at (1+`zem`)λ_Lya (where λ_Lya = 121.567 nm is the rest-frame wavelength of the Lyman-alpha transition).

## Convert y axis

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.y_convert</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>yunit</code>: Unit (of flux density)</li>
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Convert the y axis to different units.** The y axis can be expressed in different units depending on how it was calibrated (default: erg/(cm^2 s nm)). It can be converted to any unit of the same physical quantity.

## Scale y axis 🚧

## Shift to and from frame 🚧

## Rebin spectrum 🚧

## Convolve with gaussian  🚧

## Estimate resolution  🚧

## Estimate SNR 🚧
