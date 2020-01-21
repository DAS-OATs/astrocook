---
layout: default
title: Data structures
parent: Managing data
nav_order: 1
---

# Data structures
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

Astrocook manages three main data structures: *spectra*, *line lists*, and *system lists*. All three structures include a data table (actually an [Astropy Table](https://docs.astropy.org/en/stable/table/) object) and a metadata dictionary. They are formatted as FITS files and bundled into a `.acs` archive when you save a snapshot of a session. The `.acs` archive may also contain other ancillary data, as described below.

## Spectra

Each Astrocook session is based on a spectrum, i.e. a table that pairs a flux-like quantity to a wavelength-like quantity. We use the expressions “flux-like” and “wavelength-like” to be as general as possible: from the point of view of structure, it doesn't matter if the spectrum is expressed in wavelengths or frequencies, or if the flux density is calibrated in physical units or not.

The fundamental columns of a spectrum are:
- `x`: the wavelength-like independent variable;
- `xmin`, `xmax`: the interval in `x` values in which the flux-like quantity is integrated (this interval is also called a *pixel* because in some cases it maps to a physical pixel in the spectrograph detector);
- `y`: the flux-like dependent variable;
- `dy`: the error on `y`.

Other columns that frequently appear in a spectrum are:
- `conv`: a convolution of `y` with some smoothing kernel;
- `cont`: the spectral *continuum*, i.e. the component of `y` that remains after removing local features (e.g. absorption lines) and smoothing out noise;
- `resol`: the spectral resolution at `x`.
