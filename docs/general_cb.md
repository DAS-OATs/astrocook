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

### Import structure

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>astrocook.gui.GUIPanelSession.struct_import</code></td>
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

**Import a data structure from a session into the current one**. The structure to be imported is described by a string with the session number and the structure tag ( (`spec`, `lines`, `systs`) separated by a comma (e.g. `0,spec`, meaning "spectrum from session 0"). The imported structure is either replaced or appended to the corresponding one in the current session.

### Modify structures 🚧

<table>
  <tbody>
    <tr>
      <td style="vertical-align:top"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>astrocook.gui.GUIPanelSession.struct_modify</code></td>
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

**Modify a data structure using a binary operator**. An output column is computed applying a binary operator to two input columns, or an input column and a scalar. Columns are described by a string with the session number, the structure tag ( (`spec`, `lines`, `systs`), and the column name separated by a comma (e.g. `0,spec,x`, meaning "column x of spectrum from session 0").  They can be from different data structures only if they have the same length. If the output column already exists, it is overwritten.

### Apply template 🚧

### Extract region 🚧

### Convert axis 🚧

### Scale y axis 🚧

### Shift to and from frame 🚧

### Rebin spectrum 🚧

### Convolve with gaussian  🚧

### Estimate resolution  🚧

### Estimate SNR 🚧
