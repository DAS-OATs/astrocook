---
layout: default
title: Data structures
parent: Managing data
nav_order: 1
---

# Data structures
{: .no_toc}

Astrocook manages three main data structures: *spectra*, *line lists*, and *system lists*. All three structures include a data table (actually an [Astropy Table](https://docs.astropy.org/en/stable/table/) object) and a metadata dictionary. These structures and other ancillary structures are bundled into a `.acs` archive when you [save a snapshot of a session](http://localhost:4000/astrocook/gui.html#save-sessions).

---
## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

## Spectra

Each Astrocook session is based on a spectrum, i.e. a table that pairs a flux-like quantity to a wavelength-like quantity. We use the expressions “flux-like” and “wavelength-like” to be as general as possible: from the point of view of structure, it doesn't matter if the spectrum is expressed in wavelengths or frequencies, or if the flux density is calibrated in physical units or not.

[Spectrum tables](http://localhost:4000/astrocook/tables.html#spectrum-table) are saved both in FITS and ASCII format, with file names `snapshot_spec.fits/.dat` respectively (assuming `snapshot.acs` is the name of the archive containing them).

If the continuum level of the flux has been determined, a table of continuum nodes may be present, depending on how the continuum was determined (see [here](http://localhost:4000/astrocook/continuum.html#continuum-estimation) for more details). The node tables are also saved in FITS and ASCII format, with file names `snapshot_nodes.fits/.dat` respectively.


## List of lines

When you detect spectral lines (either in emission or absorption), they are formatted into a table that is added to the session.

[Line tables](http://localhost:4000/astrocook/tables.html#line-table) are saved both in FITS and ASCII format, with file names `snapshot_lines.fits/.dat` respectively.


## List of absorption systems

Absorption lines detected in quasar spectra are frequently grouped into absorption systems, i.e. sets of lines produced by different ions at the same redshift.

In our convention, a system has *one and only one* redshift. This means that e.g. different doublets at a similar redshift are treated as separated (one redshift per doublet). It also means that each component of an absorbers with a complex velocity structure (also sometimes called a "system") is considered a system in itself. [Here](series.md) you can find a complete list of the ionic transitions used to model absorption features.

[System tables](http://localhost:4000/astrocook/tables.html#system-table) are saved both in FITS and ASCII format, with file names `snapshot_systs.fits`/`.dat` respectively.

System tables only contain the parameters needed to compute the system models. The models themselves (i.e. the fitting functions, which often group together several systems into a single expression, when they need to be fitted together) are saved in a set of ASCII files `snapshot_systs_mods_NN__lines.dat`/`__group.dat`/`_left.dat`/`_right.dat`, with `NN` the ID of the model. An additional ASCII file `snapshot_systs_mods.dat` provides a list of the models. These files are handled internally via serialization by [`pickle`](https://docs.python.org/3/library/pickle.html) and [`lmfit.save_model()`/`.load_model()`](https://lmfit.github.io/lmfit-py/model.html#saving-and-loading-models) and are not meant to be accessed by the users.
