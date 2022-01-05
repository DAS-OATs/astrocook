---
layout: default
title: Tables
parent: Managing data
nav_order: 2
---

# Tables
{: .no_toc}

To display the data tables, choose `View > Spectrum table` (or `Line table` or `System table`) from the menu bar.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

## Spectrum table

The fundamental columns of a spectrum are:
- `x`: the wavelength-like independent variable;
- `xmin`, `xmax`: the interval in `x` values in which the flux-like quantity is integrated;
- `y`, `dy`: the flux-like dependent variable and its error.

The interval [`xmin`, `xmax`] is also called a *bin* (because it is the unit into which the flux information is binned) or a *pixel* (because in some cases it maps to a physical pixel in the spectrograph detector).

Other columns that frequently appear in a spectrum are:
- `y_conv`: a convolution of `y` with some smoothing kernel;
- `lines_mask`: a boolean mask of the detected lines (see [below](structures.md#list-of-lines));
- `cont`: the spectral *continuum*, i.e. the component of `y` that remains after removing local features (e.g. absorption lines) and smoothing out noise;
- `model`: a spectral *model*, i.e. a model of both the continuum and the absorption systems (see [below](structures.md#list-of-absorption-systems);
- `deabs`: a *deabsorbed* spectrum, i.e. an equivalent of `y` after the absorption lines have been removed using the model;
- `resol`: the spectral resolution at `x`;
- `fit_mask`: a boolean mask of the regions used to fit the model of the absorption systems.  

You are free to add other columns to the spectrum or edit the existing ones (as explained [here](tables.md)).


## Line table

The information to populate the table columns is directly extracted from the [spectrum](structures.md#spectra) where the lines have been detected:
- `x`: the line center, corresponding to the `x` value of a pixel in the spectrum;
- `xmin`, `xmax`: the boundaries of the interval covered by the line, also corresponding to the `x` values of two pixels in the spectrum;
- `y`, `dy`: values of `y` and `dy` at `x` from the spectrum;
- `source`: column of the spectrum that was used to detect the line.

Please note that in this case [`xmin`, `xmax`] do not map to a single pixel in the spectrum, but to a range of pixels.

Lines are not typically detected on the raw `y` column of the spectrum, because in general the noise on `y` makes it very hard to discriminate between legitimate lines and random fluctuations. This is the reason for keeping track of the `source` of the detected lines. A typical `source` value is `y_conv`, i.e. a convolution of the `y` column of the spectrum with some smoothing kernel (see [above](structures.md#spectra)) but it may be any other column.


## System table

When you detect absorption systems, they are also formatted into a table added to the session. The table contains the information to model the absorption systems:
- `func`: the function used to define the model (currently, only the Voigt function `voigt` is available);
- `series`: the list of ionic transitions that are modeled;
- `z0`: the starting redshift, typically computed from the `x` value of one or more [absorption lines](structures.md#list-of-lines);
- `z`, `dz`: the redshift of the model and its error;
- `logN`, `dlogN`: the base-10 logarithm of the column density used in the Voigt function (in cm^-2) and its error;
- `b`, `db`: the doppler parameter used in the Voigt function and its error;
- `resol`: the resolution adopted by the model;
- `chi2r`: the reduced chi-squared between the model and the data;
- `id`: the identification number of the model.

The parameters of the [Voigt profile](absorbers.md#voigt-profile-modeling) (`z`, `logN`, and `b`) may be either guess or fitted parameters. As a rule, when a system is added to the list it is also fitted to the spectrum: `z`, `logN`, and `b` are the best-fit values and `chi2r` gives an estimation of the goodness of fit.


## Browse tables 🚧

## Interact with the data 🚧
