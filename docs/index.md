# Astrocook

**Astrocook** is a Python package and GUI for analyzing quasar spectra, built on the modern Scientific Python stack (**PySide6**, **Astropy**, **Matplotlib**).  
It is designed to make spectral analysis—such as continuum fitting, absorption line identification, and Voigt profile modeling—interactive, reproducible, and easy.

## Key Features

* **Interactive GUI**: Zoom, pan, and select regions with a mouse-driven interface.  
* **Spectral Operations**: Smooth, rebin, and arithmetic operations (addition, subtraction) on spectra.  
* **Continuum Fitting**: Automatic and manual algorithms for defining quasar continua.  
* **System Identification**: Smart algorithms to identify absorption systems (Lyα, CIV, MgII, etc.).  
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

api/modules
```

## Credits

Astrocook is developed by **Guido Cupani**.  
If you use Astrocook in your research, please cite: **[Cupani et al. 2020](https://ui.adsabs.harvard.edu/abs/2020SPIE11452E..1UC/abstract)**.
