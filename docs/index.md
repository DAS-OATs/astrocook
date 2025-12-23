```{image} _static/icon_3d_HR.png
:alt: Astrocook Logo
:width: 200px
:class: no-background
:align: center
```

# Welcome to Astrocook!

**Astrocook** is a Python package and GUI for analyzing quasar spectra, built on the modern Scientific Python stack (**PySide6**, **Astropy**, **Matplotlib**).  
It is designed to make spectral analysis—such as continuum fitting, absorption line identification, and Voigt profile modeling—interactive, reproducible, and easy.

## Key Features

* **Interactive GUI**: Zoom, pan, and select regions with a mouse-driven interface.  
* **Spectral Operations**: Smooth, rebin, and arithmetic operations (addition, subtraction) on spectra.  
* **Continuum Fitting**: Automatic and manual algorithms for defining quasar continua.  
* **System Identification**: Smart algorithms to identify absorption systems (Ly-$\alpha$, CIV, MgII, etc.).  
* **Voigt Profile Fitting**: Integrated fitting engine for complex absorption features.  
* **Session Management**: Save and restore your full analysis state (undo/redo supported).

## Documentation

```{toctree}
:maxdepth: 2  
:caption: User Guide

installation
tutorials/index
```

```{toctree}
:maxdepth: 2  
:caption: API Reference

reference
```

```{toctree}
:maxdepth: 2
:caption: Credits

credits
```