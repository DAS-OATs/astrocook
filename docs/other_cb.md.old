---
layout: default
title: Other cookbook
parent: Cookbooks
nav_order: 6
---

# Other cookbook
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

###  Convert x axis
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.x_convert</code></td>
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
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "x_convert",
  "params": {
    "zem": "0",
    "xunit": "km / s"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Convert the x axis to wavelength or velocity units. The x axis can be converted to any unit of wavelength or velocity (default: nm and km/s). The conversion applies to both the spectrum and the line list. When converting to and from velocity units, the zero point is set at (1+zem)λ_Lya (where λ_Lya = 121.567 nm is the rest-frame wavelength of the Lyman-alpha transition).

###  Convert y axis
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.y_convert</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>yunit</code>: Unit of flux density</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "y_convert",
  "params": {
    "yunit": "electron / nm"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Convert the y axis to different units. The y axis can be expressed in different units depending on how it was calibrated (default: erg/(cm^2 s nm)). It can be converted to any unit of the same physical quantity. The conversion applies to both the spectrum and the line list.

###  Shift to barycentric frame
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.shift_bary</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>v</code>: Velocity in the barycentric frame (km/s)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "shift_bary",
  "params": {
    "v": "null"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Shift x axis to the barycentric frame of the solar system.

###  Shift to rest frame
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.shift_to_rf</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>z</code>: Redshift to use for shifting</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "shift_to_rf",
  "params": {
    "z": "0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Shift x axis to the rest frame.

###  Shift from rest frame
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.shift_from_rf</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>z</code>: Redshift to use for shifting</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "shift_from_rf",
  "params": {
    "z": "0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Shift x axis from rest frame to the original frame.

###  Import a data structure from a session into the current one
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.struct_import</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>struct</code>: Structure</li>
          <li><code>mode</code>: Mode (replace or append)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "struct_import",
  "params": {
    "struct": "'0,systs'",
    "mode": "'replace'"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

The structure to be imported is described by a string with the session number and the structure tag (spec, lines, systs), separated by a comma (e.g. 0,spec, meaning "spectrum from session 0"). The imported structure is either replaced or appended to the corresponding one in the current session.

###  Modify a data structure using a binary operator
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookEdit.struct_modify</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>col</code>: Output column</li>
          <li><code>expr</code>: Expression</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "struct_modify",
  "params": {
    "col": "''",
    "expr": "''"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Modify a data structure using a binary operator. An output column is computed from an expression with input columns as arguments. The expression must be parsable by AST, with columns described by a string with the session number, the structure tag (spec, lines, systs), and the column name separated by a comma (e.g. 0,spec,x, meaning "column x of spectrum from session 0"). Columns can be from different data structures only if they have the same length. If the output column already exists, it is overwritten.








