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

###  New systems from likelihood
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_new_from_like</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>series</code>: Series of transitions</li>
          <li><code>z_start</code>: Start redshift</li>
          <li><code>z_end</code>: End redshift</li>
          <li><code>dz</code>: Threshold for redshift coincidence</li>
          <li><code>modul</code>: Modulation of the error function</li>
          <li><code>thres</code>: Threshold for accepting likelihood</li>
          <li><code>distance</code>: Distance between systems in pixels</li>
          <li><code>logN</code>: Guess (logarithmic) column density</li>
          <li><code>b</code>: Guess Doppler broadening</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>chi2r_thres</code>: Reduced chi2 threshold to accept the fitted model</li>
          <li><code>dlogN_thres</code>: Column density error threshold to accept the fitted model</li>
          <li><code>refit_n</code>: Number of refit cycles</li>
          <li><code>chi2rav_thres</code>: Average chi2r variation threshold between cycles</li>
          <li><code>max_nfev</code>: Maximum number of function evaluation</li>
          <li><code>append</code>: Append systems to existing system list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_new_from_like",
  "params": {
    "series": "'Ly-a'",
    "z_start": "0",
    "z_end": "6",
    "dz": "0.0001",
    "modul": "1",
    "thres": "0.997",
    "distance": "10",
    "logN": "14",
    "b": "10",
    "resol": "null",
    "chi2r_thres": "inf",
    "dlogN_thres": "inf",
    "refit_n": "0",
    "chi2rav_thres": "0.01",
    "max_nfev": "1000",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

TBD

###  Complete systems from likelihood
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_complete_from_like</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>series</code>: Series of transitions</li>
          <li><code>series_ref</code>: Reference series of transitions</li>
          <li><code>z_start</code>: Start redshift</li>
          <li><code>z_end</code>: End redshift</li>
          <li><code>binz</code>: Bin size to group existing redshifts</li>
          <li><code>dz</code>: Threshold for redshift coincidence</li>
          <li><code>modul</code>: Modulation of the error function</li>
          <li><code>thres</code>: Threshold for accepting</li>
          <li><code>distance</code>: Distance between systems in pixels</li>
          <li><code>logN</code>: Guess (logarithmic) column density</li>
          <li><code>b</code>: Guess Doppler broadening</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>chi2r_thres</code>: Reduced chi2 threshold to accept the fitted model</li>
          <li><code>dlogN_thres</code>: Column density error threshold to accept the fitted model</li>
          <li><code>refit_n</code>: Number of refit cycles</li>
          <li><code>chi2rav_thres</code>: Average chi2r variation threshold between cycles</li>
          <li><code>max_nfev</code>: Maximum number of function evaluation</li>
          <li><code>append</code>: Append systems to existing system list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_complete_from_like",
  "params": {
    "series": "'all'",
    "series_ref": "null",
    "z_start": "0",
    "z_end": "6",
    "binz": "0.01",
    "dz": "0.0001",
    "modul": "1",
    "thres": "0.997",
    "distance": "10",
    "logN": "14",
    "b": "10",
    "resol": "null",
    "chi2r_thres": "inf",
    "dlogN_thres": "inf",
    "refit_n": "0",
    "chi2rav_thres": "0.01",
    "max_nfev": "1000",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

TBD

###  New systems from line list
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_new_from_lines</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>series</code>: Series of transitions</li>
          <li><code>z_start</code>: Start redshift</li>
          <li><code>z_end</code>: End redshift</li>
          <li><code>dz</code>: Threshold for redshift coincidence</li>
          <li><code>logN</code>: Guess (logarithmic) column density</li>
          <li><code>b</code>: Guess Doppler broadening</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>chi2r_thres</code>: Reduced chi2 threshold to accept the fitted model</li>
          <li><code>dlogN_thres</code>: Column density error threshold to accept the fitted model</li>
          <li><code>refit_n</code>: Number of refit cycles</li>
          <li><code>chi2rav_thres</code>: Average chi2r variation threshold between cycles</li>
          <li><code>max_nfev</code>: Maximum number of function evaluation</li>
          <li><code>append</code>: Append systems to existing system list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_new_from_lines",
  "params": {
    "series": "'Ly-a'",
    "z_start": "0",
    "z_end": "6",
    "dz": "0.0001",
    "logN": "14",
    "b": "10",
    "resol": "null",
    "chi2r_thres": "inf",
    "dlogN_thres": "inf",
    "refit_n": "0",
    "chi2rav_thres": "0.01",
    "max_nfev": "1000",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Add and fit Voigt models to a line list, given a redshift range.

###  Complete systems
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_complete</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>series</code>: Series of transitions</li>
          <li><code>dz</code>: Threshold for redshift coincidence</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>avoid_systs</code>: Avoid adding transitions over systems already fitted</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_complete",
  "params": {
    "series": "'all'",
    "dz": "0.0001",
    "resol": "null",
    "avoid_systs": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Add candidate transitions to fitted systems.

###  Find candidate systems
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.cands_find</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>series</code>: Series of transitions</li>
          <li><code>z_start</code>: Start redshift</li>
          <li><code>z_end</code>: End redshift</li>
          <li><code>dz</code>: Threshold for redshift coincidence</li>
          <li><code>resol</code>: Resolution</li>
          <li><code>avoid_systs</code>: Avoid finding candidates over systems already detected</li>
          <li><code>append</code>: Append systems to existing system list</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "cands_find",
  "params": {
    "series": "'all'",
    "z_start": "0",
    "z_end": "6",
    "dz": "0.0001",
    "resol": "null",
    "avoid_systs": "true",
    "append": "true"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Cross-match line wavelengths with known transitions to find candidate systems.

###  Improve systems
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_improve</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>impr_n</code>: Number of improve cycles</li>
          <li><code>refit_n</code>: Number of refit cycles</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_improve",
  "params": {
    "impr_n": "3",
    "refit_n": "0"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Improve systems adding components to reduce residuals

###  Fit systems
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_fit</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>refit_n</code>: Number of refit cycles</li>
          <li><code>chi2rav_thres</code>: Average chi2r variation threshold between cycles</li>
          <li><code>max_nfev</code>: Maximum number of function evaluation</li>
          <li><code>sel_fit</code>: Selective fit (only new systems will be fitted)</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_fit",
  "params": {
    "refit_n": "3",
    "chi2rav_thres": "0.01",
    "max_nfev": "1000",
    "sel_fit": "false"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Fit all Voigt model from a list of systems.

###  Clean system list
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_clean</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>chi2r_thres</code>: Reduced chi2 threshold to accept the fitted model</li>
          <li><code>dlogN_thres</code>: Column density error threshold to accept the fitted model</li>
          <li><code>max_nfev</code>: Maximum number of function evaluation</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_clean",
  "params": {
    "chi2r_thres": "2.0",
    "dlogN_thres": "1.0",
    "max_nfev": "1000"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Clean systems from a list by rejecting systems with reduced chi2 and/or error on column density above a given threshold

###  Recreate the models
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.mods_recreate</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        –
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "mods_recreate",
  "params": {
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Recreate the models from the current system list.

###  Estimate SNR of systems
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_snr</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        –
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_snr",
  "params": {
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Estimate the signal-to-noise ratio of systems as the median flux/flux error ratio in the group interval.

###  Select systems
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_select</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>series</code>: Series</li>
          <li><code>z_min</code>: Minimum redshift</li>
          <li><code>z_max</code>: Maximum redshift</li>
          <li><code>logN_min</code>: Minimum (logarithmic) column density</li>
          <li><code>logN_max</code>: Maximum (logarithmic) column density</li>
          <li><code>b_min</code>: Minimum Doppler broadening</li>
          <li><code>b_max</code>: Maximum Doppler broadening</li>
          <li><code>col</code>: Other column</li>
          <li><code>col_min</code>: Minimum of other column</li>
          <li><code>col_max</code>: Maximum of other column</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_select",
  "params": {
    "series": "'any'",
    "z_min": "0.0",
    "z_max": "10.0",
    "logN_min": "10.0",
    "logN_max": "22.0",
    "b_min": "1.0",
    "b_max": "100.0",
    "col": "null",
    "col_min": "null",
    "col_max": "null"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Select systems based on their Voigt and fit parameters. A logical `and` is applied to all conditions.

###  Extract systems based on components
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.comp_extract</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>num</code>: Number of components</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "comp_extract",
  "params": {
    "num": "1"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Extract systems with less than a given number of components

###  Merge a system into the current system
        
<table>
  <tbody>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>
      <td style="vertical-align:top"><code>CookbookAbsorbers.systs_merge</code></td>
    </tr>
    <tr>
      <td style="vertical-align:top"><strong>Parameters</strong></td>
      <td style="vertical-align:top">
        <ul>
          <li><code>num1</code>: Row of the current system</li>
          <li><code>num2</code>: Row of the system to be merged</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>
      <td style="vertical-align:top"><pre>
{
  "cookbook": "cb",
  "recipe": "systs_merge",
  "params": {
    "to_row": "0",
    "from_rows": "[1]"
  }
}    </pre></td>
    </tr>
  </tbody>
</table>

Merged systems appear as a single entry in the compressed system table.

