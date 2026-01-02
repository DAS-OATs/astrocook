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

## 2. Navigating the Plot

The central plot is interactive. You can explore the details of your spectrum using the toolbar located directly above the graph (standard Matplotlib controls):

- **Pan**: Click the "Cross arrow" icon to click and drag the spectrum.
- **Zoom**: Click the "Magnifying glass" icon to draw a rectangle and zoom into a specific feature.
- **Home**: Click the "Home" icon to reset the view to the full range.

:::{tip}
Hover your mouse over any data point or absorption line. A tooltip will appear displaying the precise Wavelength, Flux, and information about identified components or regions!
:::

(getting-started-customizing)=
## 3. Customizing the View (Right Sidebar)

Look to the right side of the window. This is the **Plot Controls** panel. Here you can tweak how the data is displayed without altering the data itself.

### Toggles

You can show or hide specific layers of the plot:

- **1-sigma error**: Hides the grey error shading.
- **Normalization**: Hides the continuum level (dashed black line).
- **Absorption model**: Hides the Voigt profile fits (solid red line).
- **Systems**: Hides the vertical tick marks indicating identified systems.
- **Aux. Column**: Select an additional data column (e.g., `abs_mask`, `running_std`) to overlay on the plot.
- **Strong emission lines**: Shows markers for major emission lines (e.g., Ly-$\alpha$, CIV). 
  :::{note}
  This requires defining the emission redshift first (see [Session Properties](#editing-properties)).
  :::
- **Show Iso-Mag Grid**: Overlays curves of constant AB magnitude. 
  :::{note} 
  This requires flux calibration (see [Flux Operations](flux_operations.md)).
  :::

### Axis Controls

- **Units**: You can switch the X-axis display between `nm`, `Angstrom`, and `micron` on the fly.
- **Normalize F**: Plots the flux relative to the continuum (showing $f = F / F_\mathrm{cont}$).
- **Show SNR**: Plots the Signal-to-Noise Ratio instead of flux.
    * **Error Col.**: When showing SNR, you can select which column to use as the noise estimate (default is `dF`).
- **Logarithmic $\lambda$**: Switches the X-axis to a logarithmic scale.

### Axis Limits

You can manually set the precise bounds of the plot view. Enter values for **$\lambda$ Min/Max** and **F Min/Max** and press Enter to apply them. These boxes automatically update when you zoom or pan interactively.

### The Redshift Cursor

This is a handy tool for quick visual inspection.

1. In the **Redshift Cursor** section, enter a transition name (e.g., `Ly_a`, `CIV`, `MgII`) and a redshift `z`.
2. Check the **Show Cursor Lines** box.
3. Vertical dashed lines will appear on the plot indicating where that transition should be found. You can update the `z` value to slide the cursor along the spectrum.

```{image} ../_static/getting_started_plot_controls.png
:alt: The Plot Controls sidebar
:align: center
```

**Interactive Setting:** You can **Right-click** directly on the plot background to open a context menu. Select **Set emission redshift to z=...** to instantly update the session's $z_{em}$ based on the cursor's current position.

:::{note}
If the **Zoom** or **Pan** tools are active in the toolbar, the right-click shortcut is overridden. In this case, hold `Ctrl` (or `Cmd` on Mac) while right-clicking to access the context menu.
:::


## 4. Session Management (Left Sidebar)

The panel on the left is your **Session List**. Each file you load becomes a "Session".

### Switching Sessions

If you loaded multiple spectra, click on their names in this list to switch the main view instantly.

### Session Info & Properties

Right-click on any session in the list to open the context menu. Select **View Info** to open the Session Inspector. Here you can see (and edit) critical metadata:

- **Object Name**
- **Emission Redshift ($z_\mathrm{em}$)**: Crucial for many recipes like the Ly-$\alpha$ forest analysis.
- **Resolution ($R$)**: Essential for accurate Voigt profile fitting.

### The Log Scripter

Astrocook automatically records every operation you perform (smoothing, fitting, editing). Right-click a session and select **View Log** to open the **Log Scripter**.

This dialog allows you to:
* **Review History:** See a readable list of every recipe run on this session.
* **Undo/Redo:** Step back through your analysis history.
  :::{tip}
  You can access Undo/Redo also from the **Edit** menu (see [Editing Data](editing_data.md)) or using the standard key combinations `Ctrl+Z`/`Shift+Ctrl+Z` (`Cmd+Z`/`Shift+Cmd+Z` on Mac).
  :::
* **Run All:** Re-run the entire analysis pipeline from the original raw data.
* **Save Script:** Export your workflow as a Python-compatible JSON script for reproducibility.

### Saving and Closing

When you have finished your analysis (or if you want to save your progress):

1. **Save**: Go to **File > Save Session...** or right-click the session in the list. This saves your work as an `.acs2` file, preserving all your continuum fits and line lists.
2. **Close**: You can close a specific session via the right-click menu, or close the active session via **File > Close Session**.

```{image} ../_static/getting_started_session_list.png
:alt: The Session List with context menu
:align: center
```

## What's Next?

Now that you are comfortable with the interface, you are ready to manipulate your data. Check out the [Editing Data](editing_data.md) tutorial to learn how to clean up your spectra.