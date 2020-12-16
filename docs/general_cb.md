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

## Open session

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>GUIPanelSession._on_open</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>path</code>: Path to an <code>.acs</code>, <code>.fits</code>, or <code>.json</code> file</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><code>
      {
        "cookbook": "_panel_sess",
        "recipe": "_on_open",
        "params": {
          "path": "XXX"
        }
      },
    </code></td>
    </tr>
  </tbody>
</table>

**Open a new session.** Session can be opened from previously saved Astrocook archive (`.acs`), a FITS file (`.fits`), or a JSON file (`.json`). In the last case, the session is built by running the workflow in the JSON file, and the workflow is inherited in the session log.

## Equalize sessions

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>GUIPanelSession.equalize</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>xmin</code>: Minimum wavelength (nm)</li>
          <li><code>xmax</code>: Maximum wavelength (nm)</li>
          <li><code>_sel</code>: Selected sessions</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><code>
    {
      "cookbook": "_panel_sess",
      "recipe": "equalize",
      "params": {
        "xmin": "XXX",
        "xmax": "XXX",
        "_sel": [
          XXX,
          XXX
        ]
      }
    }
    </code></td>
    </tr>
  </tbody>
</table>

**Equalize the flux level of one session to another one.** The last-selected session is equalized to the first-selected one. The equalization factor is the ratio of the median flux within the specified wavelength interval. When the recipe is called from the GUI, the `_sel` parameter is defined automatically by clicking on the two input sessions.

## Combine sessions

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>GUIPanelSession.combine</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>name</code>: Name of the output session</li>
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Combine two or more sessions.** A new session is created, with a new spectrum containing all entries from the spectra of the combined sessions. Other objects from the sessions (line lists, etc.) are discarded. When the recipe is called from the GUI, the `_sel` parameter is defined automatically by clicking on the input sessions.

## Import structure

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
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

**Import a data structure from a session into the current one.** The structure to be imported is described by a string with the session number and the structure tag (`spec`, `lines`, `systs`), separated by a comma (e.g. `0,spec`, meaning "spectrum from session 0"). The imported structure is either replaced or appended to the corresponding one in the current session.

## Modify structures

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
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
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
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

**Extract a spectral region.** The region between a minimum and a maximum wavelength is extracted from the data structures in the current session (these include the selected spectral range with all the lines, and the absorption systems that fall within). A new session with the extracted data structures is created.

## Convert x axis

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
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

**Convert the x axis to wavelength or velocity units.** The x axis can be converted to any unit of wavelength or velocity (default: nm and km/s). The conversion applies to both the spectrum and the line list. When converting to and from velocity units, the zero point is set at (1+`zem`)λ_Lya (where λ_Lya = 121.567 nm is the rest-frame wavelength of the Lyman-alpha transition).

## Convert y axis

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
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

**Convert the y axis to different units.** The y axis can be expressed in different units depending on how it was calibrated (default: erg/(cm^2 s nm)). It can be converted to any unit of the same physical quantity. The conversion applies to both the spectrum and the line list.

## Scale y axis

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookGeneral.y_scale</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>fact</code>: Multiplicative factor</li>
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Scale the y axis by a constant factor.** The spectrum and the line list are rescaled in place, without starting a new session.

## Scale y axis to median 🚧

## Shift to and from frame 🚧

## Rebin spectrum

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
        </ul>
      </td>
    </tr>
  </tbody>
</table>

**Rebin a spectrum with a given step.** The step can be expressed in any unit of wavelength or velocity. Start and end wavelength may be specified, e.g. to align the rebinned spectrum to other spectra. If start or end wavelength are `None`, rebinning is performed from the first to the last wavelength of the input spectrum. A new session is created with the rebinned spectrum. Other objects from the old session (line lists, etc.) are discarded.

## Convolve with gaussian  🚧

## Estimate resolution  🚧

## Estimate SNR 🚧

## Refresh the GUI

_refresh(self, init_cursor=False, init_tab=True, autolim=True,
                 autosort=True, _xlim=None):
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>GUI._refresh</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>init_cursor</code>: Initialize system cursor</li>
          <li><code>init_tab</code>: Initialize tables</li>
          <li><code>autolim</code>: Automatically set limits to the plot axes</li>
          <li><code>autosort</code>: Automatically sort tables</li>
          <li><code>_xlim</code>: Limits for the plot x axis</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><code>
      {
        "cookbook": "",
        "recipe": "_refresh",
        "params": {
          "autosort": false
        }
      }
    </code></td>
    </tr>
  </tbody>
</table>

**Refresh the GUI.** This recipe is designed for internal use. The user should call it only at the end of a workflow, to update the visualization of the data in the GUI. `autosort` should be set to `false` not to interfere with other possible sorting of the tables in the workflow. Other parameters should be disregarded.
