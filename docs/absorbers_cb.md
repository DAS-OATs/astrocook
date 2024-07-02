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
      <td style="vertical-align:top"><code>CookbookAbsorbers.lines</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>kind</code>: Kind of lines (`abs` or `em`)</li>
          <li><code>prominence</code>: Prominence of lines (as in `scipy.signal.find_peaks`)</li>
          <li><code>append</code>: Append lines to existing line list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "lines",
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
