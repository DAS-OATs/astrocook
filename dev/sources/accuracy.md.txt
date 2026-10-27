# Accuracy

How close are the results of Astrocook to the truth? This page reports
what the automatic chain gives on synthetic spectra with known
parameters. The numbers refer to version 2.2.0. They are a guide to
where the results can be trusted and where they should be checked, not
a guarantee for real data.

## How it is measured

The mock spectra are made by `astrocook.mocks`, which does not use the
fitting code of Astrocook. Each spectrum covers 104-160 nm in the rest
frame of the quasar, from the Lyman-beta forest to beyond the CIV
emission. It has:

- a power-law continuum with emission lines;
- a Lyman-alpha forest, with the whole Lyman series for logN >= 15;
- metal systems of CIV, SiIV and MgII on the red side, with their H I;
- the Lyman limit, and a noise that follows the flux;
- in two spectra out of five, also a DLA and a sub-DLA in the forest.

The grid has 135 spectra: S/N 10, 30 and 100 per pixel; resolving power
R = 8900, 45000 and 140000 (2.5 pixels per resolution element); emission
redshift 2.5, 3.5 and 4.5; five random realisations of each.

On each spectrum the chain runs with the default parameters of the
dialogs:

1. **Continuum > Auto-estimate Continuum**;
2. **Absorbers > Identify Lines**;
3. **Absorbers > Auto-detect Doublets**.

A fitted line is matched to a true one within 20 km/s. The script is
`benchmarks/validate_mocks.py` in the source repository; the full grid
takes about 25 minutes on four cores.

## Continuum

The table gives the ratio of the estimated to the true continuum, as
the median over the spectra. In brackets is the 5th percentile, which
shows how low the continuum gets in the worst regions.

| | Ly-alpha forest | Red side |
|---|---|---|
| All | 0.943 (0.830) | 1.002 (0.996) |
| S/N 10 / 30 / 100 | 0.873 / 0.938 / 0.980 | 1.006 / 1.002 / 1.001 |
| z_em 2.5 / 3.5 / 4.5 | 0.971 / 0.938 / 0.889 | 1.002 |

On the red side the continuum is accurate. In the forest it is
systematically low, because the unabsorbed pixels are few. The bias is
larger at low S/N and at high redshift. Check the continuum in the
forest, and correct it by hand where needed (see
[Continuum Fitting](tutorials/continuum_fitting.md)).

## Identification of metal doublets

The table gives the purity (the fraction of the candidates that are
real) and the completeness. The completeness counts only the systems
whose weaker line is detected at 3 sigma or more.

| | Purity | Completeness |
|---|---|---|
| CIV | 0.99 | 0.87 |
| SiIV | 0.99 | 0.69 |
| MgII | 0.99 | 0.90 |

The completeness depends on the S/N:

| S/N | CIV | SiIV | MgII |
|---|---|---|---|
| 10 | 0.76 | 0.53 | 0.79 |
| 30 | 0.89 | 0.71 | 0.92 |
| 100 | 0.94 | 0.81 | 0.99 |

It also depends on the significance of the weaker line of the doublet
(its equivalent width over its error):

| Significance | 3-5 | 5-10 | > 10 |
|---|---|---|---|
| Completeness | 0.10 | 0.44 | 0.91 |

A candidate that is found is almost always real. Weak systems, however,
are mostly missed: below a significance of about 10, look for them by
hand.

## Column densities

### Total column of a system

The total column of a system is the sum of the fitted components within
its velocity range, compared with the true sum. It is reliable:

| R | Fitted systems | Median difference | Spread |
|---|---|---|---|
| 8900 | 247 of 360 | -0.044 dex | 0.19 dex |
| 45000 | 291 of 360 | -0.002 dex | 0.075 dex |
| 140000 | 302 of 360 | -0.002 dex | 0.031 dex |
| All | 840 of 1080 | -0.006 dex | 0.077 dex |

The spread is half the interval between the 16th and the 84th
percentile.

### Single components

For an isolated component, the parameters and their errors are right.
Its pulls are the difference from the truth divided by the error given
by the fit. Their width is 1.10 for z, 1.03 for logN and 1.15 for b,
close to the ideal value of 1.

The decomposition of a system with several close components is not
reliable: the automatic fit often uses fewer components than the true
ones, and each fitted component then stands for several. The
parameters of these components, and their errors, should not be read
as those of single clouds: use the total column of the system, or
refine the fit by hand in the
[System Inspector](tutorials/absorption_systems.md).

## Fit quality

Of the components fitted by the pipeline, 11 out of 1181 have no true
counterpart. 22 are flagged by the fit, for a parameter at a bound or
for uncertainties that could not be computed (see *Components with
warnings of the fit* in [Absorption Systems](tutorials/absorption_systems.md)).
The median reduced chi2 is 1.25. It is higher at R = 8900 (1.6), where
the components are less resolved.
