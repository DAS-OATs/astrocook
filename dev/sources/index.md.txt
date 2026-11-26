```{image} _static/icon_3d_HR.png
:alt: Astrocook Logo
:width: 200px
:class: no-background
:align: center
```

# Welcome to Astrocook!

*A thousand recipes to cook a spectrum.*


**Astrocook** is a Python package and GUI for analyzing quasar spectra, built on the modern Scientific Python stack (**PySide6**, **Astropy**, **Matplotlib**).  
It is designed to make spectral analysis—such as continuum fitting, absorption line identification, and Voigt profile modeling—interactive, reproducible, and easy.

The tagline explains the name: every operation of Astrocook (smoothing, estimating the continuum, identifying the systems, fitting the profiles, ...) is a **recipe**, a step with its own parameters, available from the menus, from Python scripts and from AI agents, and recorded in the history of the session. An analysis is a sequence of recipes, applied to a spectrum.

## Key Features

* **Interactive GUI**: Zoom, pan, and select regions with a mouse-driven interface.  
* **Spectral Operations**: Smooth, rebin, and arithmetic operations (addition, subtraction) on spectra.  
* **Continuum Fitting**: Automatic and manual algorithms for defining quasar continua.  
* **System Identification**: Smart algorithms to identify absorption systems (Ly-$\alpha$, CIV, MgII, etc.).  
* **Voigt Profile Fitting**: Integrated fitting engine for complex absorption features, with the instrumental resolution of each pixel and warnings on the components whose result is not reliable.  
* **Session Management**: Save and restore your full analysis state (undo/redo supported).

## Documentation

```{toctree}
:maxdepth: 1  
:caption: User Guide

installation
settings
tutorials/index
science
reference
```

```{toctree}
:maxdepth: 1
:caption: Info

credits
changelog
```