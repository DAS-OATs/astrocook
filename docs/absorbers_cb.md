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
