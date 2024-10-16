---
layout: default
title: Synthetic cookbook
nav_order: 4
---

# Synthetic cookbook
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Synthetic spectrum from structures

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookSynthetic.spec_from_struct</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>x</code>: Expression for wavelength-like array</li>
          <li><code>y</code>: Expression for flux-like array</li>
          <li><code>dy</code>: Expression for the error on y</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "spec_from_struct",
  "params": {
    "x": "'0,spec,x'",
    "y": "'0,spec,y'",
    "dy": "'0,spec,y'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Create a synthetic spectrum from existing structures (a wavelenght-like array and a flux-like array). The structure expressions must be parsable by AST, with columns described by a string with the session number, the structure tag (spec, lines, systs), and the column name separated by a comma (e.g. 0,spec,x, meaning "column x of spectrum from session 0"). A gaussian noise is added to the spectrum to match a given signal-to-noise ratio. A new session is created with the synthetic spectrum.

###  Synthetic spectrum from systems

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookSynthetic.spec_from_systs</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>x</code>: Expression for wavelength-like array</li>
          <li><code>dy</code>: Expression for the error on y</li>
          <li><code>sess</code>: Number of the session with the systems</li>
          <li><code>resol</code>: Resolution</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "spec_from_systs",
  "params": {
    "x": "'0,spec,x'",
    "dy": "'0,spec,dy'",
    "sess": "'0'",
    "resol": "null"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Create a synthetic spectrum from a list of systems taken from an existing session.

###  Synthetic spectrum from random systems

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookSynthetic.spec_from_systs_random</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>n</code>: Number of systems</li>
          <li><code>series</code>: Series of transitions</li>
          <li><code>z_min</code>: Minimum redshift</li>
          <li><code>z_max</code>: Maximum redshift</li>
          <li><code>z_seed</code>: Seed for random sampling in [z_min, z_max]</li>
          <li><code>logN_min</code>: Minimum (logarithmic) column density</li>
          <li><code>logN_max</code>: Maximum (logarithmic) column density</li>
          <li><code>logN_seed</code>: Seed for random sampling in [logN_min, logN_max]</li>
          <li><code>b_min</code>: Minimum Doppler broadening</li>
          <li><code>b_max</code>: Maximum Doppler broadening</li>
          <li><code>b_seed</code>: Seed for random sampling in [b_min, b_max]</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>snr</code>: Signal-to-noise ratio</li>
          <li><code>append</code>: Append systems to existing system list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "spec_from_systs_random",
  "params": {
    "n": "n",
    "series": "'Ly-a'",
    "z_min": "0",
    "z_max": "6",
    "z_seed": "null",
    "logN_min": "10",
    "logN_max": "18",
    "logN_seed": "null",
    "b_min": "1.0",
    "b_max": "100.0",
    "b_seed": "null",
    "resol": "null",
    "snr": "null",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Create a synthetic spectrum from a list of systems with random redshifts, column density, and Doppler broadening.
