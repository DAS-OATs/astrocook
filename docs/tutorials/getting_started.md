# Getting Started

Welcome to Astrocook V2! This guide will walk you through the basics of the interface: loading data, exploring your spectra, and managing your analysis sessions.

## 1. Loading a Spectrum

Let's begin by loading some data. Astrocook supports various formats, including standard FITS files and its own archive formats (`.acs`, `.acs2`).

1.  Launch Astrocook. You will be greeted by the welcome screen.
2.  Go to the menu bar and select **File > Open Session...** (or press `Ctrl+O` / `Cmd+O`).
3.  Select your spectrum file. You can select multiple files at once if you want to load several observations simultaneously.

Once loaded, the interface will come alive. The central area displays your spectrum, while two sidebars provide control over your data and visualization.

```{image} ../_static/getting_started_main_view.png
:alt: The main Astrocook interface with a spectrum loaded
:align: center
```

(getting-started-navigating)=
## 2. Navigating the Plot

The central plot is interactive. You can explore the details of your spectrum using the toolbar below the graph. The tool that is on (Pan, Zoom or Select) is highlighted in blue; Back and Forward are grey when there is no view to go to. The keys in brackets act after a click on the plot:

- **Home** (`H`): reset the view to the full range.
- **Back** and **Forward** (`←`, `→`): the previous and next views.
- **Pan** (`P`): drag the spectrum; drag with the right button to zoom.
- **Zoom** (`Z`): draw a rectangle to zoom into a specific feature.
- **Select** (`X`): the "Range" icon (a double arrow between two bars); drag along the wavelength axis to select a region. When you release the mouse, a menu offers what to do with it: **Split region** (a new session with the region only), **Mask region...**, **Rebin region...** (each opens its recipe dialog with the region filled in, so you can check it before running) and **Zoom to region**. Press `Esc` or click elsewhere to discard the region. The tool stays on until you click its icon again or choose Pan or Zoom; it is not available when the x axis is in velocity.
- **Save**: save the figure to an image file.

`Esc` switches off the tool that is on. The position of the mouse (wavelength, flux and rest-frame wavelength) is shown at the right end of the status bar, at the bottom of the window; its left end says how the last recipe went (e.g. *Auto-fit Ly-alpha forest: done in 12.3 s, 44 → 977 components*, and how many warnings it gave). **View > Show Toolbar** hides the toolbar (and shows it again), to give the plot more room; its keys still work, the position of the mouse stays in the status bar, and the choice is remembered.

:::{tip}
Hover your mouse over any data point or absorption line. A tooltip will appear displaying the precise Wavelength, Flux, and information about identified components or regions!
:::

```{image} ../_static/getting_started_region_select.png
:alt: A region of the spectrum selected with the Select tool, and the menu of the actions on it
:align: center
```

(data-inspector)=
### The Data Inspector

For precise pixel-level analysis, you can open the **Data Inspector**.

1.  **Right-click** on any point in the spectrum plot.
2.  Select **Inspect Data (near ...)** from the context menu.

TThis opens a floating window displaying the raw values centered on the pixel you clicked. The table headers use scientific notation for clarity:
* **λ**: Central wavelength of the pixel.
* **λmin / λmax**: The lower and upper boundaries of the pixel bin.
* **F**: Flux density.
* **dF**: 1-sigma error estimate.
* **Other columns**, containing the continuum and additional information, when present (such as the resolution column `resol`, mentioned below).

The table supports infinite scrolling, automatically loading more data as you move up or down.

(getting-started-customizing)=
## 3. Customizing the View (Right Sidebar)

Click the cog icon at the top right of the plot to open the **Plot Controls** panel (click it again to close it). Here you can tweak how the data is displayed without altering the data itself.

The panel floats over the plot, so that opening it does not change the aspect of the spectrum. To cover as little as possible, it is only as tall as its contents, and it fades out while the mouse is on the plot: you can see the effect of a change through it. It comes back as soon as the mouse returns to it. Its sections (**Plot Elements**, **View Controls**, **Axis Limits**, **Redshift Cursor**, **Edit Continuum**) open and close with a click on their title, and Astrocook remembers which ones you left open.

### Plot Elements

You can show or hide specific layers of the plot. Next to each one, a short sample of the line shows how it is drawn: these samples are the legend of the plot, which has none of its own. The auxiliary column chosen in **Auxiliary column** has its sample inside the menu, and each other session shown (see below) next to its name in the list of the sessions.

- **1-sigma error**: Shows or hides the error of the flux (grey).
- **Show error as line**: Checked by default: the error array ($dy$) is drawn as a line of its own, useful for inspecting the noise structure in high-resolution spectra. Unchecked, the error is drawn as a band around the flux.
- **Normalization**: Hides the continuum level (dashed black line).
- **Absorption model**: Hides the Voigt profile fits (solid green line).
- **Systems**: Hides the tick marks of the fitted components at the bottom of the plot, in two rows: HI (red) below, metals (grey) above.
- **Auxiliary column**: Select an additional data column (e.g., `abs_mask`, `running_std`) to overlay on the plot.
- **Strong emission lines only**: Shows only the major emission lines (e.g., Lyα, C IV), named in spectroscopic notation; unchecked, all the lines of the atomic data, fainter. 
  :::{note}
  This requires defining the emission redshift first (see [Setting Session Properties](#editing-properties)).
  :::
- **Show iso-magnitude grid**: Overlays curves of constant AB magnitude. 
  :::{note} 
  This requires flux calibration (see [Flux Calibration](#flux-calibration)).
  :::

### Axis Controls

- **Units**: You can switch the X-axis display between `nm`, `Angstrom`, and `micron` on the fly.
- **Normalize F**: Plots the flux relative to the continuum (showing $f = F / F_\mathrm{cont}$) (this control becomes active only when a `cont` column is present).
- **Show S/N**: Plots the signal-to-noise ratio instead of flux.
    * **Error column**: When showing the S/N, you can select which column to use as the noise estimate (default is `dF`).
- **Logarithmic $\lambda$**: Switches the X-axis to a logarithmic scale.
- **Logarithmic F**: Switches the Y-axis to a logarithmic scale, e.g. to see the continuum and the emission lines over several orders of magnitude, or the **iso-magnitude grid** as straight lines. Pixels with $F \le 0$ (saturated absorption, noise) are drawn at the bottom of the plot. In the normalized view the axis starts from 0.01, so saturated lines are shown down to 1% of the continuum. A non-positive **F Min** is not applied while the scale is logarithmic.

### Axis Limits

You can manually set the precise bounds of the plot view. Enter values for **$\lambda$ Min/Max** and **F Min/Max** and press Enter to apply them. These boxes automatically update when you zoom or pan interactively.

### The Redshift Cursor

This is a handy tool for quick visual inspection.

1. In the **Redshift Cursor** section, enter a transition name (e.g., `Ly_a`, `CIV`, `MgII`) and a redshift `z`.
2. Check the **Show cursor** box.
3. Vertical dashed lines will appear on the plot indicating where that transition should be found, with their redshift written at the top of the first one. The cursor follows the mouse over the plot; you can also type the `z` value to place it.

```{image} ../_static/getting_started_plot_controls.png
:alt: The Plot Controls sidebar
:align: center
```

**Interactive Setting:** You can **Right-click** directly on the plot background to open a context menu. Select **Set Emission Redshift to z=...** to instantly update the session's $z_{em}$ based on the cursor's current position.

:::{note}
If the **Zoom** or **Pan** tools are active in the toolbar, the right-click shortcut is overridden. In this case, hold `Ctrl` (or `Cmd` on Mac) while right-clicking to access the context menu.
:::

### Dark Mode

By default Astrocook follows the light or dark appearance of your system (on macOS, *System Settings > Appearance*), and switches with it while it is open. **View > Appearance** chooses *Light* or *Dark* for Astrocook alone, or *System* to follow the system again; the choice is remembered. The plots stay white in both, like the figures of a paper, so that the colours of the lines are always the same.

```{image} ../_static/getting_started_dark_mode.png
:alt: The main window in dark mode, with the plot white
:align: center
```

### The Continuum Editor

Located at the bottom of the sidebar, this panel allows you to interactively refine the continuum shape using spline knots. See the [Interactive Manual Editing](#continuum-manual) in the Continuum Fitting tutorial for a detailed guide on using the slider, spacing box, and editing tools.

## 4. Session Management (Left Sidebar)

Click the list icon at the top left of the plot to open the **Session List**. Each file you load becomes a "Session". Like the Plot Controls, the list floats over the plot, is as tall as the sessions it holds and fades out while the mouse is on the plot. Long names are shortened in the middle; hover over a name to read it in full.

### Switching Sessions

If you loaded multiple spectra, click on their names in this list to switch the main view instantly.

### Visualizing Multiple Spectra (Ghost Overlay)

You can easily compare multiple spectra simultaneously directly in the main plot window:

* **How to use:** Hold `Ctrl` (or `Cmd` on macOS) and click on multiple sessions in the left-hand Session List. You can also use `Shift+Click` to select a continuous range.
* **What you'll see:** Your primary active session remains fully featured (showing models, continuums, and system markers). The additional selected sessions will appear underneath as "ghost" spectra, plotted in distinct, colorblind-friendly colors with semi-transparent flux and error bands. A sample of its colour appears next to the name of each of them in the session list, to help you identify them.
* **Smart Scaling:** The ghost spectra are fully integrated with your View Controls! If you toggle **Normalize F**, all visible spectra will automatically rescale to match the active view.

```{image} ../_static/getting_started_ghost_spectra.png
:alt: Overlay of multiple ghost spectra in the main plot window
:align: center
```

### Session Info & Properties

Right-click on any session in the list to open the context menu. This menu offers several powerful tools:

- **View Info**: Opens the *Session Info* dialog. Here you can see (and edit) critical metadata:
  - **Object**
  - **Emission Redshift ($z_\mathrm{em}$)**: Crucial for many recipes like the Ly-$\alpha$ forest analysis.
  - **Resolution ($R$)**: Essential for accurate Voigt profile fitting.
  :::{note}
  Setting a global resolution here will automatically generate a pixel-by-pixel `resol` column across your entire spectrum, which you can verify in the Data Inspector.
  :::

- **Duplicate**: Creates a complete copy of the current session. This allows you to branch your analysis (e.g., to test two different continuum fits on the same data).

### The Session Log

Astrocook automatically records every operation you perform (smoothing, fitting, editing). Right-click a session and select **View Log** to open the **Session Log**.

This dialog allows you to:
* **Review History:** See a readable list of every recipe run on this session.
* **Undo/Redo:** Step back through your analysis history.
  :::{tip}
  You can access Undo/Redo also from the **Edit** menu (see [Editing Data](editing_data.md)) or using the standard key combinations `Ctrl+Z`/`Ctrl+Y` (`Cmd+Z`/`Cmd+Shift+Z` on Mac).
  :::
* **Run All:** Re-run the entire analysis pipeline from the original raw data.
* **Save Script:** Export your workflow as a Python-compatible JSON script for reproducibility.

### View Menu Tools

The **View** menu also provides access to specialized analysis windows:

* **View Session Log**: Opens the Session Log (see above).
* **View System Inspector**: Opens the advanced interface for analyzing absorption lines. *See the [Absorption Systems](#system-inspector) tutorial for a full guide.*
* **View Identifications**: Opens a summary table of all candidate lines found by the `identify_lines` recipe.
* **View Script Console**: Opens a Python console for running custom scripts on the current session data.

### Saving and Closing

When you have finished your analysis (or if you want to save your progress):

1. **Save**: Go to **File > Save Session...** or right-click the session in the list. This saves your work as an `.acs2` file, preserving all your continuum fits, line lists and the session log (the operations shown in the Log Scripter).
2. **Close**: You can close a specific session via the right-click menu, or close the active session via **File > Close Session**.

```{image} ../_static/getting_started_session_list.png
:alt: The Session List with context menu
:align: center
```

## 5. Exporting and Importing Data

Astrocook allows you to easily exchange data with external tools (like spreadsheet managers, Topcat, or custom Python scripts) using standard CSV files.

### Export to ASCII
To save your current data in a human-readable format, go to **File > Export to ASCII...**. 
A dialog will appear allowing you to select exactly what to export:
* **Spectrum (spec)**: Saves the Wavelength ($x$), Flux ($y$), Error ($dy$), Continuum ($cont$), and any other active spectral columns into a `_spec.csv` file.
* **System List (systems)**: Saves your identified absorption lines (including Redshift, Column Density, Doppler Parameter, etc.) into a `_systems.csv` file.

You can choose to export entire structures or pick individual columns. The suggested filename will automatically adjust based on your selection, and Astrocook will warn you if you are about to overwrite existing files.

### Import from ASCII
If you have manipulated your exported `_spec.csv` data externally, you can bring it back into Astrocook by going to **File > Import from ASCII...**.
* Astrocook will match the column headers in your CSV to the columns in the active session and overwrite them.
* **Safety first:** The imported file must have the exact same number of rows (pixels) as your current spectrum. Because this operates through the recipe pipeline, the import is fully undoable via the Log Scripter or `Ctrl+Z` (`Cmd+Z`).

## What's Next?

Now that you are comfortable with the interface, you are ready to manipulate your data. Check out the [Editing Data](editing_data.md) tutorial to learn how to clean up your spectra.