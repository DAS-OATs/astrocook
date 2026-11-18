# Identification as a comparison between hypotheses

Design document of phase 7, PR 7c.2 (`ROADMAP.md`, phase 7, block 7c).
Written on 18 November 2026 and read by Guido the same day: his
decisions, and the proposals added after his questions (accepted the
same day), are in the last section. No production code comes with this document: the experiment of
section 4 is `benchmarks/identification_hypotheses.py`.

The goal of the block (phase 7): *every absorbed region accounted for,
identified or explicitly unidentified, with alternatives and a
credibility, tellurics recognised*. The proposal, in one paragraph: the
absorbed regions are found by a detector of their own (not by the clipping
of the continuum), grouped by the transitions that can link them; for each
group a short list of alternative hypotheses (metal X at z, Ly-a, an
unidentified line, telluric, BAL or DLA from 7b), each a forward model of
all its lines convolved with the resolution, measured on the same pixels by
`evaluate_model` (PR 7a.5); a cheap evaluation for all (N and b from the
apparent optical depth), a short fit for the finalists; the credibility of
each the weight exp(-ΔIC/2) of an information criterion, BIC or the SpIC of
Webb et al. (2021) by the experiment of section 4; the alternatives kept
with their weights in the data structure of 7c.3.

## 1. The current algorithm and its numbers

`identify_lines` (`astrocook/core/spectrum_lines.py`, recipe
`absorbers.identify_lines`, default multiplets `CIV,SiIV,MgII,Ly_ab`):

1. *Regions.* The absorbed pixels are the mask `abs_mask` of
   `find_absorbed` (`find_absorbed_regions`, `core/spectrum_operations.py`):
   the pixels below the running mean of the iterative clipping of the
   continuum by more than kappa = 2 sigma (window 5000 km/s in the forest,
   400 on the red side). Contiguous pixels (at least 3) are a region;
   regions closer than `merge_dv` = 10 km/s are merged.
2. *Kinematic candidates.* For each doublet of `STANDARD_MULTIPLETS` with
   exactly two lines (CIV, SiIV, MgII, ...; the longer multiplets are
   never searched, F93 item 2), two regions are a candidate when they
   overlap once one is shifted by the lag of the doublet; the metals on the
   red side of the Ly-a emission only (in the forest the threshold is
   1.1, impossible), the `Ly_ab` pairs in the forest only; no candidate
   beyond z_em by more than 3000 km/s (F55).
3. *Score.* `rate_doublet_candidate`: the R^2 of F_1 = F_2^gamma (gamma =
   f_1 lambda_1 / f_2 lambda_2) on the normalised flux in the overlap
   window (padded by 30 km/s), in the logarithmic frame of the lags since
   PR 7c.1; a candidate is accepted above 0.5. The relation holds for the
   optical depth before the convolution: after it, and for unresolved or
   saturated lines, it holds only approximately.
4. *Ly-a by default.* Every region of the forest without an accepted
   candidate is labelled Ly-a with score 0 (F51), whatever the multiplets
   asked; the dropped candidates are listed with their score in
   `details.dropped_candidates` (F119, PR 7-ter of 6-ter).
5. *De-blending.* Each identified region gives one candidate system per
   peak of 1 - F (prominence max(0.02, 3 sigma)), with N from the AOD around
   the peak and b from its width; `doublets_auto` and `components_auto`
   then fit them.

The pairs are not stored (each region keeps its [transition, score]
pairs, F93 item 1); there are no alternatives, no credibility, no
telluric hypothesis, and the "unidentified" label does not exist: a region
of the red side without a candidate is simply not identified.

**Its numbers** (PRs 7a.2 and 7c.1, reference
`benchmarks/results/validate_mocks_robustness_baseline.json`):

| | `full` | `robustness` |
|---|---|---|
| CIV purity / completeness | 0.993 / 0.870 | 0.859 / 0.809 |
| SiIV purity / completeness | 0.987 / 0.857 | 0.823 / 0.779 |
| MgII purity / completeness | 0.988 / 0.913 | 0.864 / 0.808 |
| CIV completeness by sigma of the weaker line, 3-5 / 5-10 / > 10 | 0.12 / 0.55 / 0.96 | |
| SiIV, same | 0.17 / 0.45 / 0.93 | |
| MgII, same | - / 0.29 / 0.96 | |
| False candidates on tellurics (all species) | | 172 (MgII 57 of its 64 false) |
| False candidates on associated absorbers | | 49 CIV, 39 SiIV (7a.2) |

The pixels (7a.2): of the detectable absorbed pixels, 73-88% carry only
the Ly-a by default, about 1% an identification that was tested, 10-22% are
in no region at all (78% of the BAL pixels, 96% of the telluric ones). The
Ly-a/Ly-b pairs (`validate_lyab_pairs.py`, 7c.1): purity 0.84, completeness
0.087 over the components with a detectable Ly-b (0.27 at z_em 2.5, 0.05 at
4.5): the dense forest is merged into long regions, and a pair is found per
pair of regions, not per component.

What the numbers say: the identification of the metals on the red side is
pure and, above 10 sigma, complete; it loses most systems between 3 and 10
sigma (F57), where the experiment of section 4 shows that the information
is there (a doublet at 3-5 sigma of its weaker line is told from noise in
every trial, with the true continuum); the forest is not identified but
assumed; nothing of the BALs and tellurics is recognised.

## 2. The literature, and what to reuse

References from memory are marked (*); their details are to be checked
before they enter the docs or the paper.

**Automatic Voigt fitters.**
- *AutoVP* (Davé, Hernquist, Weinberg & Katz 1997, ApJ 477, 21): the
  forest of simulations decomposed automatically, components added at the
  minima of the residuals until a chi2 criterion is met. All lines are
  Ly-a: no identification.
- *VPFIT* (Carswell & Webb 2014, ascl:1408.015) with *RDGEN*: the
  interactive standard; the identification is the user's, VPFIT adds and
  removes components by chi2 and the significance of each. The `fort.26`
  format is the import of 7a.3.
- *GVPFIT* (Bainbridge & Webb 2017, MNRAS 468, 1639)(*) and *AI-VPFIT*
  (Lee, Webb, Carswell & Milaković 2021, MNRAS 504, 1787)(*): a genetic
  algorithm, then a non-deterministic placement of components, with the
  number of components decided by an information criterion; Webb et al.
  (2021, MNRAS 501, 2268) introduce the SpIC for that decision and find
  AICc overfitting and BIC underfitting (section 5 of their paper). Their
  setting is one species set known in advance: they choose the components,
  not the identification. What to reuse: the SpIC (in `evaluate_model`
  since 7a.5) and the idea that a criterion, not a threshold on a score,
  decides.
- *VIPER* (Gaikwad et al. 2017, MNRAS)(*): automatic Voigt fitting of the
  forest with the number of components chosen by AICc. *BayesVP* (Liang &
  Kravtsov 2017, arXiv:1710.09852)(*): MCMC fits with the Bayesian
  evidence for the number of components (the "dynesty for important
  cases" of section 5 of the roadmap).

**Doublet searches in surveys.** The CIV and MgII catalogues of SDSS
(e.g. Nestor, Turnshek & Rao 2005, ApJ 628, 637, MgII; Cooksey et al.
2013, ApJ 763, 37, CIV; Seyffert et al. 2013, ApJ 779, 161, MgII; Zhu &
Ménard 2013, ApJ 770, 130, MgII)(*) search the doublets on the spectrum
normalised by a smooth continuum (a spline, or a PCA/NMF model of the
quasar), with a filter matched to the line profile (a Gaussian of the
instrumental width), and accept a pair on two thresholds (typically 4-5
sigma on the stronger line, 2.5-3 on the weaker), the doublet ratio
within its physical range (1 to 2 for thin to saturated) and the
separation within a tolerance; completeness and purity measured by
injecting simulated doublets into real spectra. What to reuse: the
detection on the normalised flux with a matched filter (section 5), the
two thresholds of a pair as a cheap pre-filter of the candidates, the
injection of doublets into real spectra as a measure (the twins of 7a.3
are the same idea).

**Classifiers.** DLAs: Gaussian processes on the forest (Garnett et al.
2017, MNRAS 472, 1850)(*), convolutional networks (Parks et al. 2018,
MNRAS 476, 1151)(*); MgII and CIV: networks trained on the SDSS catalogues.
They are trained on low-resolution survey spectra, and give a probability
without a model: what to reuse is the habit of reporting a probability,
not the classifiers (the DLAs belong to 7b, the metals of high-resolution
spectra have too few labelled examples).

**Tellurics.** Masks of the bands (the simplest; e.g. the masks of the
radial-velocity pipelines), templates scaled by airmass, and models fitted
on the spectrum: *Molecfit* (Smette et al. 2015, A&A 576, A77; Kausch et
al. 2015, A&A 576, A78), *TelFit* (Gullikson et al. 2014, AJ 148, 53)(*),
on the LBLRTM radiative transfer; SkyCalc (Noll et al. 2012, A&A 543,
A92; Jones et al. 2013, A&A 560, A91) for the model of 7a.1. The fit of the
telluric transmission is 7b.6; for the identification what matters is
that the telluric hypothesis is a model like the others (section 7).

**What no tool does** (as far as known): the identification of the
metals of a high-resolution spectrum as a comparison of forward models of
alternative identifications with a credibility. The closest are the
interactive habits of the observers (look for the other lines of a
candidate ion; prefer a doublet to two forest lines when the ratio is
right) and the automatic search of the SDSS catalogues, which the proposal
combines.

## 3. The hypotheses

### 3.1 Groups

The absorbed regions (section 5) are the nodes of a graph; two regions
are linked when a transition pair of the requested species can put them at
the same redshift (within the tolerance of the kinematic search, which
exists: `_find_kinematic_doublet_candidates`), the two lines of a doublet,
two lines of a longer multiplet (section 8), Ly-a and Ly-b (in the forest),
and the lines of different ions of one system at the same z (CIV and SiIV,
MgII and FeII; section 3.5). The connected components of the
graph are the *groups* (the "blocks" of section 5 of the roadmap): the
alternatives of a group do not change the evaluation of another, so the
groups are evaluated one by one, in parallel if needed. These groups are
not the groups of the fits (`SystemList.get_connected_group`, the
components whose profiles overlap): an identification group joins regions
that one hypothesis can explain together, also far apart in wavelength;
the components of the accepted hypotheses are then fitted in the groups of
the fits, as today (Guido, 18 November 2026).

In the dense forest the graph percolates: every region is linked to some
other by the Ly-a/Ly-b lag or by a metal lag. Two limits keep the groups
small: (a) the edges are kept only for candidates that pass a cheap
pre-filter, the R^2 score of today with a low threshold (decided: Guido,
18 November 2026; the threshold to measure in 7c.4), (b) a group
above a size (to agree, e.g. 6 regions) is split at its weakest edges,
those of the lowest pre-filter significance. The regions of the forest
that no edge reaches are groups of one region, with the hypotheses Ly-a
(when asked, section 6) and unidentified.

### 3.2 Alternatives

For each group a short list of hypotheses, each explaining *every* region
of the group (otherwise the chi2 punishes the one that leaves a region
unmodelled, not the one that is wrong):
- *metal X at z*: the components of the species (as many as the peaks of
  the de-blending), all the lines of X in the coverage (section 8);
- *Ly-a*: one component per peak, with its Ly-b where covered (the higher
  orders not in the model of the identification: section 12); in the forest only, or beyond it within the margin of F55;
- *unidentified*: one free line per peak, modelled as a Voigt profile of a
  placeholder transition (the Ly-a profile at its wavelength, k = 3), with
  no consequence elsewhere; on the red side it is the alternative of a
  metal; in the forest it has the model of Ly-a without its Ly-b;
- *telluric*: the telluric transmission alone (section 7);
- *BAL* or *DLA*: the region inside a trough or a damped profile of the
  emission model of 7b; these regions are taken out of the groups and
  labelled by 7b (the common decision on the vocabulary, last section).

A hypothesis mixes labels (in the alias scenario of section 4: a CIV at
z' and a Ly-a for the line it leaves). The enumeration per group is: each
accepted kinematic candidate of the group as the anchor of one hypothesis,
the other regions completed by the default label of their side (Ly-a in
the forest when asked, unidentified otherwise), plus the hypothesis "all
default". A group of n candidates gives n + 1 hypotheses, not 2^n; the
pairs of candidates that do not share a region can be combined in one
hypothesis (two metal systems in one group) when both win alone.

### 3.3 Evaluation and credibility

Each hypothesis is a system list (the components of the hypothesis plus
the fixed rest of the spectrum), measured by `evaluate_model` on the
pixels of the group: the union of the windows around every line of every
hypothesis of the group (the same pixels for all, as the comparison
requires). The model is convolved with the resolution, so the relation
F_1 = F_2^gamma is replaced by the forward model where it does not hold
(unresolved lines, saturation), and a multiplet of any length is a model
like a doublet. Every component of a hypothesis is counted
(`min_influence = 0`): a hypothesis pays for all its parameters.

Two stages:
1. *Without a fit*, for every hypothesis: N from the apparent optical depth
   (AOD) of the lines (the largest over the lines of a multiplet, which
   corrects the saturation of the stronger), b from the second moment of
   N(v) with the instrumental width subtracted, z from the candidate. A
   screen: the hypotheses far behind (ΔIC above a threshold, e.g. 20) are
   dropped.
2. *A short fit* (300 evaluations, the bounds of the pipelines) of the
   finalists from their AOD start, and the comparison on the fitted ones.
   The experiment (section 4) shows why the second stage is needed: the
   AOD start of a doublet is a poor model of saturated lines, and the
   comparison without a fit favours the hypothesis with more free lines
   (the doublet against two Ly-a: 76-94% right by AOD, 97-99% after the
   fit, SpIC_H); and that it is cheap (a few ms per fit of an isolated
   system, section 10).

The credibility of hypothesis i of a group is its weight
w_i = exp(-ΔIC_i / 2) / Σ_j exp(-ΔIC_j / 2), ΔIC_i = IC_i - min_j IC_j
(Akaike weights for AIC, the approximate posterior for BIC; for the SpIC
the same form, by analogy, as Webb et al. 2021 do not define weights). The
criterion: BIC or SpIC_H by the experiment (section 4); proposed SpIC_H as
the default, BIC reported next to it, and the choice measured again in
7c.4 on the full mocks (decided: SpIC_H, BIC reported; Guido, 18 November 2026).

The weights assume a chi2 that is calibrated: errors that are right and
independent pixels. The reductions are rebinned, so neighbouring pixels
are correlated, and the error spectrum can be off by some per cent; a
chi2 inflated by a factor s inflates every ΔIC by s, and the weights
become overconfident (a ΔIC of 10, w = 0.99, is 5 with s = 2, w = 0.92).
The mocks of the experiment of section 4 have independent pixels, so its
weights are the optimistic case. Proposed: the weights tempered by the
chi2 per degree of freedom of the accepted model on the pixels of the
group (ΔIC / max(1, chi2_nu)), or by an effective number of pixels from
the correlation of the noise when it is known; which one, measured by the
calibration of F164 (section 11) on the twins and the development
spectra, where the correlation is real (after the conversation of PR #798,
18 November 2026).

### 3.4 A prior on the incidence of the ions

The weights above assume equal prior odds of the hypotheses. In the
forest they are not equal: a line of logN 13 is far more likely a Ly-a
than a CIV. The prior odds of two hypotheses enter as
ΔIC → ΔIC - 2 ln(P_1 / P_2).

The right odds are not the ratio of the line densities dN/dz of the ions:
the Ly-a hypothesis of a doublet candidate needs *two* Ly-a lines at the
lag of the doublet, the second within the tolerance of the search. At
z = 3 the forest of the mocks has dN/dz = 1700 for 12 < logN < 17
(`astrocook/mocks/populations.py`, set to give the mean transmission of
Faucher-Giguère et al. 2008): about 0.023 lines per km/s, so a second line
within ±5 km/s by chance has a probability of about 0.2 (0.07 above logN
13); CIV systems above logN 12.5 have dN/dz of order 5-10 at z ~ 2-3
(D'Odorico et al. 2010; Cooksey et al. 2013 for the stronger ones)(*).
The odds of the Ly-a pair against the CIV are then of order 1700 × 0.2 / 10
≈ 30 for weak lines, of order 1-10 above logN 13: not negligible, and
dependent on z, logN and the tolerance. A proper prior would integrate
over the column densities (the evidence, not the maximum likelihood),
which the criteria do not.

The experiment gives the effect (section 4.3): with odds 10 for the Ly-a,
SpIC still keeps 82% of the true CIV doublets in the forest, with odds 100
only 16% (BIC 98% and 86%). Proposed:
- no prior in the prototype (odds 1), the option in the recipe with its
  sources (a table of dN/dz per ion and z range in the package, cited);
- the metals in the forest stay off by default (as today: a prior of
  zero), and open on request; where a metal of the forest has another
  ion of the same system on the red side (CIV in the forest with MgII
  outside; SiIV with CIV), the link of section 3.1 makes the group, and
  the other ion is the evidence the prior cannot give;
- the prior of the mocks must not come from the mocks: the forest of
  `lya_forest` has its own dN/dz (1700 (1 + z)^2.5 / 4^2.5), and a prior
  taken from it would measure the generator, not the method.

### 3.5 The sweep of the anchored systems: rarer and single lines

(Added after Guido's reading, 18 November 2026.) The anchors (doublets
and the pairs of section 8) find the systems; most other lines of a
system are too weak, too rare or single (CII 1334, SiIII 1206, OI 1302,
AlII 1670, the weaker FeII and SiII lines) to be found by a pair. They are
identified as consequences of the anchored systems, in a second pass:
- for each accepted system at z, and for each species of the table with
  lines in the coverage (all species, not only the requested anchors), a
  hypothesis "Y at z": components at the redshifts of the anchor's
  components (the kinematics of ions of different ionisation differ: z
  free within the width of the anchor's profile, b free), N from the AOD
  at the predicted positions;
- evaluated on the pixels of its predicted lines against the labels those
  pixels already have (unidentified, the assumed Ly-a, another
  hypothesis): accepted when it lowers the criterion of the groups
  involved, i.e. when it explains regions left unexplained at the cost of
  its parameters; the z of a known system is the prior that makes a single
  line credible, where alone it would be unidentified;
- where nothing is there, the same evaluation gives the upper limit of
  the column of Y (F156, 7d.7), from the same pixels;
- an accepted species becomes an edge of the graph: the regions it
  explains join the group of its system, so a later comparison sees them
  together (a "Ly-a" of the forest that is the SiII 1260 of a system of
  the red side).
Cost: some tens of systems by some tens of species with a line in the
coverage, each a cheap evaluation on a few windows: of the order of the
first pass. The same sweep, with Ly-a and Ly-b, tests the H I of a
metal system (its lines in the forest).

## 4. The experiment: controlled pairs of alternatives

`benchmarks/identification_hypotheses.py`. Each trial is a mock (flat
continuum, the true continuum given, a window of a few nm, no forest) with
a known truth and two hypotheses, one true, measured by `evaluate_model`
on the same pixels (the union of ±200 km/s around every line of either),
without a fit (AOD) and after a short fit. The scenarios:

| scenario | truth | true hypothesis (k) | false hypothesis (k) |
|---|---|---|---|
| `doublet_vs_lya` | CIV doublet, logN 11.8-14.8, b 4-20 | CIV (3) | two Ly-a on its lines (6) |
| `doublet_vs_lyab` | the same | CIV (3) | two Ly-a with their Ly-b (6), the Ly-b pixels added |
| `coincidence` | two Ly-a at the CIV lag (±5 km/s), logN 12.3-14.5 and b 15-40 each, independent | two Ly-a (6) | CIV (3) |
| `alias` | CIV doublet | CIV (3) | CIV with its 1548 on the true 1550, and a Ly-a on the true 1548 (6) |
| `detection` | CIV doublet (the low columns below the detection limit) | CIV (3) | no absorber (0) |

The grid: S/N 10, 30, 100 by R 8900, 45000, 100 trials per scenario and
cell (3000 trials); the significance of each trial is that of its weaker
line (CIV 1550, or the weaker Ly-a of a coincidence), with the rule of
`validate_mocks.is_detectable`. A trial picks the truth when the criterion
of the false hypothesis is larger.

**Limits of the design.** The lines are isolated: no blends with the
forest, the true continuum (no continuum error), the Ly-b of the
`doublet_vs_lyab` falls on a clean continuum (in a spectrum, in the Ly-b
forest: this is its best case). The `coincidence` draws its two Ly-a
independently, so some pairs have the ratio of a doublet by chance (then
the data cannot tell, and a prior is the only help). The criteria are
measured on the maximum likelihood, not on the evidence.

### 4.1 Results

3000 trials in 4.5 minutes on 4 cores (`--trials 100`, the default 40).
Fraction of the trials that pick the truth, after the short fit, by the
significance of the weaker line:

| scenario | sigma | n | chi2 | aic | aicc | bic | spic | spic_a | spic_b |
|---|---|---|---|---|---|---|---|---|---|
| `doublet_vs_lya` | <3 | 152 | 0.16 | 0.88 | 0.91 | 0.99 | 0.97 | 0.98 | 0.83 |
| `doublet_vs_lya` | 3-5 | 43 | 0.23 | 0.91 | 0.91 | 1.00 | 0.98 | 1.00 | 0.93 |
| `doublet_vs_lya` | 5-10 | 87 | 0.08 | 0.94 | 0.95 | 0.99 | 0.97 | 0.98 | 0.98 |
| `doublet_vs_lya` | >10 | 318 | 0.04 | 0.92 | 0.93 | 1.00 | 0.99 | 0.92 | 1.00 |
| `doublet_vs_lyab` | <3 | 201 | 0.30 | 0.91 | 0.91 | 0.99 | 0.97 | 0.97 | 0.90 |
| `doublet_vs_lyab` | 3-5 | 46 | 0.46 | 0.89 | 0.91 | 1.00 | 1.00 | 1.00 | 1.00 |
| `doublet_vs_lyab` | 5-10 | 85 | 0.86 | 0.98 | 0.98 | 1.00 | 0.99 | 0.99 | 1.00 |
| `doublet_vs_lyab` | >10 | 268 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| `coincidence` | <3 | 80 | 0.95 | 0.82 | 0.82 | 0.78 | 0.80 | 0.74 | 0.84 |
| `coincidence` | 3-5 | 45 | 0.96 | 0.89 | 0.82 | 0.73 | 0.82 | 0.82 | 0.84 |
| `coincidence` | 5-10 | 101 | 0.99 | 0.91 | 0.89 | 0.86 | 0.87 | 0.88 | 0.87 |
| `coincidence` | >10 | 374 | 1.00 | 0.98 | 0.98 | 0.95 | 0.97 | 0.98 | 0.97 |
| `alias` | <3 | 175 | 0.30 | 0.92 | 0.92 | 0.99 | 0.97 | 0.99 | 0.89 |
| `alias` | 3-5 | 50 | 0.52 | 0.96 | 0.98 | 1.00 | 1.00 | 1.00 | 1.00 |
| `alias` | 5-10 | 68 | 0.93 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| `alias` | >10 | 307 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| `detection` | <3 | 146 | 0.92 | 0.60 | 0.60 | 0.36 | 0.52 | 0.47 | 0.60 |
| `detection` | 3-5 | 68 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| `detection` | 5-10 | 68 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| `detection` | >10 | 318 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

The same without a fit (N and b from the AOD), four criteria:

| scenario | sigma | n | chi2 | aic | bic | spic |
|---|---|---|---|---|---|---|
| `doublet_vs_lya` | <3 | 152 | 0.37 | 0.76 | 0.91 | 0.94 |
| `doublet_vs_lya` | 3-5 | 43 | 0.60 | 0.86 | 0.93 | 0.86 |
| `doublet_vs_lya` | 5-10 | 87 | 0.61 | 0.74 | 0.84 | 0.76 |
| `doublet_vs_lya` | >10 | 318 | 0.72 | 0.76 | 0.79 | 0.76 |
| `doublet_vs_lyab` | <3 | 201 | 0.41 | 0.85 | 0.95 | 0.92 |
| `doublet_vs_lyab` | 3-5 | 46 | 0.72 | 0.91 | 1.00 | 0.96 |
| `doublet_vs_lyab` | 5-10 | 85 | 0.75 | 0.87 | 0.95 | 0.87 |
| `doublet_vs_lyab` | >10 | 268 | 0.95 | 0.96 | 0.97 | 0.96 |
| `coincidence` | <3 | 80 | 0.97 | 0.84 | 0.81 | 0.81 |
| `coincidence` | 3-5 | 45 | 0.93 | 0.84 | 0.78 | 0.82 |
| `coincidence` | 5-10 | 101 | 0.98 | 0.93 | 0.92 | 0.93 |
| `coincidence` | >10 | 374 | 0.97 | 0.95 | 0.94 | 0.95 |
| `alias` | <3 | 175 | 0.57 | 0.85 | 0.98 | 0.95 |
| `alias` | 3-5 | 50 | 0.82 | 0.96 | 1.00 | 0.96 |
| `alias` | 5-10 | 68 | 0.88 | 0.96 | 0.97 | 0.96 |
| `alias` | >10 | 307 | 0.99 | 0.99 | 0.99 | 0.99 |
| `detection` | <3 | 146 | 0.71 | 0.38 | 0.21 | 0.32 |
| `detection` | 3-5 | 68 | 1.00 | 1.00 | 0.96 | 0.99 |
| `detection` | 5-10 | 68 | 1.00 | 1.00 | 1.00 | 1.00 |
| `detection` | >10 | 318 | 1.00 | 1.00 | 1.00 | 1.00 |

After the fit, the trials above 3 sigma by cell:

| scenario | S/N | R | n | aic | bic | spic |
|---|---|---|---|---|---|---|
| `doublet_vs_lya` | 10 | 8900 | 47 | 0.96 | 1.00 | 0.98 |
| `doublet_vs_lya` | 10 | 45000 | 69 | 0.91 | 1.00 | 0.97 |
| `doublet_vs_lya` | 30 | 8900 | 70 | 0.97 | 1.00 | 1.00 |
| `doublet_vs_lya` | 30 | 45000 | 83 | 0.90 | 1.00 | 0.99 |
| `doublet_vs_lya` | 100 | 8900 | 83 | 0.93 | 0.99 | 0.99 |
| `doublet_vs_lya` | 100 | 45000 | 96 | 0.90 | 1.00 | 0.97 |
| `doublet_vs_lyab` | 10 | 8900 | 38 | 1.00 | 1.00 | 1.00 |
| `doublet_vs_lyab` | 10 | 45000 | 55 | 0.98 | 1.00 | 1.00 |
| `doublet_vs_lyab` | 30 | 8900 | 57 | 0.98 | 1.00 | 1.00 |
| `doublet_vs_lyab` | 30 | 45000 | 75 | 0.96 | 1.00 | 0.99 |
| `doublet_vs_lyab` | 100 | 8900 | 77 | 1.00 | 1.00 | 1.00 |
| `doublet_vs_lyab` | 100 | 45000 | 97 | 0.98 | 1.00 | 1.00 |
| `coincidence` | 10 | 8900 | 50 | 0.76 | 0.66 | 0.70 |
| `coincidence` | 10 | 45000 | 81 | 0.95 | 0.84 | 0.89 |
| `coincidence` | 30 | 8900 | 89 | 0.94 | 0.88 | 0.90 |
| `coincidence` | 30 | 45000 | 100 | 0.99 | 0.98 | 0.99 |
| `coincidence` | 100 | 8900 | 100 | 1.00 | 1.00 | 1.00 |
| `coincidence` | 100 | 45000 | 100 | 1.00 | 1.00 | 1.00 |
| `alias` | 10 | 8900 | 45 | 1.00 | 1.00 | 1.00 |
| `alias` | 10 | 45000 | 61 | 1.00 | 1.00 | 1.00 |
| `alias` | 30 | 8900 | 61 | 1.00 | 1.00 | 1.00 |
| `alias` | 30 | 45000 | 80 | 0.99 | 1.00 | 1.00 |
| `alias` | 100 | 8900 | 83 | 0.99 | 1.00 | 1.00 |
| `alias` | 100 | 45000 | 95 | 1.00 | 1.00 | 1.00 |
| `detection` | 10 | 8900 | 44 | 1.00 | 1.00 | 1.00 |
| `detection` | 10 | 45000 | 68 | 1.00 | 1.00 | 1.00 |
| `detection` | 30 | 8900 | 67 | 1.00 | 1.00 | 1.00 |
| `detection` | 30 | 45000 | 88 | 1.00 | 1.00 | 1.00 |
| `detection` | 100 | 8900 | 88 | 1.00 | 1.00 | 1.00 |
| `detection` | 100 | 45000 | 99 | 1.00 | 1.00 | 1.00 |

The median weight of the truth (fit, above 3 sigma) is 1.00 for every
criterion in every scenario except `doublet_vs_lya`: chi2 0.27, AIC 0.88,
AICc 0.91, BIC 1.00, SpIC_H 0.97, SpIC_A 0.91, SpIC_B 0.99.

### 4.2 What the results say

- *chi2 alone cannot identify.* Two Ly-a lines with free N and b contain
  the doublet as a special case, so their chi2 is lower in 96% of the
  trials above 10 sigma (`doublet_vs_lya`): any choice by the fit alone
  prefers the hypothesis with more lines. The penalty decides.
- *The nested case follows the theory.* For `doublet_vs_lya` the chi2
  gained by the 3 extra parameters is distributed about as chi2 with 3
  degrees of freedom: AIC (penalty 6) keeps the doublet in 1 - P(chi2_3 >
  6) ≈ 0.89 of the trials (measured 0.88-0.94), BIC (3 ln n ≈ 16 for 200
  pixels) in nearly all (0.99-1.00). SpIC_H is between them (0.97-0.99):
  its R_a is smaller than n for weak lines, so its penalty is too.
- *The coincidences are the hard case, and BIC is weakest there.* Two
  unrelated Ly-a lines at the lag of the doublet are taken for a CIV in
  9-27% of the trials between 3 and 10 sigma (BIC 0.73-0.86, SpIC_H
  0.82-0.87, AIC 0.89-0.91), at S/N 10 and R 8900 in 24-34%. The misses are
  the pairs whose columns happen to differ by about the doublet ratio:
  for SpIC_H the median |ΔlogN - log 2| of the misses is 0.22 dex, of the
  hits 0.89; the data cannot tell those apart, and only a prior (or the
  Ly-b, or another ion) can. The penalty that protects the doublet in
  `doublet_vs_lya` is what favours it here: the choice of the criterion
  is a trade between the two, and SpIC_H is the balanced one (as Webb et
  al. 2021 find for the number of components).
- *The forward model of the whole series is the strongest evidence.* With
  the Ly-b of the interlopers in the model (`doublet_vs_lyab`, its best
  case: a clean Ly-b), BIC and SpIC keep the doublet in 99-100% of the trials
  above 3 sigma (AIC 89-100%), even chi2 above 10 sigma; the Ly-b of the Ly-a
  hypothesis predicts an absorption that is not there.
- *The alias* (the wrong pairing of a kinematic search) is rejected in
  89-99% of the trials below 3 sigma by the penalised criteria, and in
  all above 5 sigma.
- *Detection is not limited by the comparison.* Between 3 and 5 sigma of
  its weaker line the doublet beats "no absorber" in every trial by every
  criterion, with the true continuum; below 3 sigma the penalised criteria
  say "no absorber" in 40-64% of the trials, as they should. The 0.12-0.17
  completeness of the benchmark at 3-5 sigma (F57) is lost in the
  detection of the regions, before any comparison (section 5).
- *The evaluation without a fit is a screen, not a judge.* The AOD
  parameters of a doublet (N the largest of its lines, one b) are a poor
  model of saturated lines: in `doublet_vs_lya` above 10 sigma SpIC_H
  picks the doublet in 76% of the trials by AOD and 99% after the fit; in
  `alias` and `detection` the AOD is nearly as good as the fit. The AOD
  can drop the alternatives far behind, the finalists need the fit.
- *Resolution and S/N.* Above 3 sigma the cells differ little, except the
  coincidences, which become harder at low S/N and low R (fewer pixels
  across a line, less shape to tell two lines from a doublet).

Proposed from these numbers: SpIC_H as the criterion, BIC reported next to
it; the Ly-b modelled with the Ly-a wherever it is covered; the short fit
of the finalists. To be measured again in 7c.4 on the full mocks, with the
forest, blends and the estimated continuum, which this experiment leaves
out.

### 4.3 The prior

Fraction of the trials (fit, above 3 sigma) that pick the truth when the
prior odds of the Ly-a hypothesis against the CIV one are added as
2 ln(odds) to its criterion:

| scenario | odds Ly-a : CIV | AIC | BIC | SpIC_H |
|---|---|---|---|---|
| `doublet_vs_lya` | 1 | 0.92 | 1.00 | 0.98 |
| `doublet_vs_lya` | 10 | 0.36 | 0.98 | 0.82 |
| `doublet_vs_lya` | 100 | 0.01 | 0.86 | 0.16 |
| `doublet_vs_lyab` | 1 | 0.98 | 1.00 | 1.00 |
| `doublet_vs_lyab` | 10 | 0.93 | 1.00 | 0.97 |
| `doublet_vs_lyab` | 100 | 0.83 | 0.99 | 0.88 |
| `coincidence` | 1 | 0.96 | 0.92 | 0.93 |
| `coincidence` | 10 | 0.99 | 0.93 | 0.97 |
| `coincidence` | 100 | 0.99 | 0.97 | 0.99 |
| `alias` | 1 | 1.00 | 1.00 | 1.00 |
| `alias` | 10 | 0.96 | 1.00 | 0.98 |
| `alias` | 100 | 0.92 | 1.00 | 0.92 |

The prior moves the trade of section 4.2: odds of 10 recover most of the
coincidences (SpIC_H 0.93 → 0.97) at the cost of 16% of the true doublets
when only the Ly-a is modelled, 3% when its Ly-b is; odds of 100 lose most
of the doublets of the forest under SpIC_H unless the Ly-b is modelled.
With the odds of section 3.4 (order 1-30, by z and logN) the prior matters
in the forest, and the Ly-b, or another ion of the system, matters more.

## 5. Detection decoupled from the continuum, and the scoring

**Detection** (F57). Today the regions are the pixels that the clipping
of the continuum rejects: a weak line pulls the running mean down by a
fraction of its depth, and its pixels do not reach 2 sigma; the
completeness at 3-5 sigma is 0.12-0.17. The experiment shows that the
comparison itself is not the limit: a CIV doublet whose weaker line is at
3-5 sigma is preferred to "no absorber" in every trial by every criterion
(section 4.1, `detection`). Proposed:
- the detection on the flux normalised by the final continuum (of 7b, or
  the current one: last section), with a bank of matched filters: for
  each pixel the significance of a Gaussian absorption of FWHM the
  instrumental one combined with b = 5, 15, 40 km/s (three widths), S =
  Σ w_i (1 - F_i) / σ_i^2 / sqrt(Σ w_i^2 / σ_i^2), w the filter; a line
  where S > 3 for any width, its extent where S > 1 (the region);
- the doublets (and multiplets) searched also on the joint significance
  of their lines at the lag, the filter of the weaker scaled by the ratio
  of f lambda: a doublet whose lines are at 3 and 1.5 sigma has a joint
  significance of about 3.4, and is tested; the look-elsewhere effect (the
  number of independent positions and lags) is what the threshold of 3
  must face, and the purity of the benchmark measures it;
- the regions of the detection replace `abs_mask` in the identification
  (`abs_mask` stays the continuum's mask for the continuum).

**Scoring** (F119). The R^2 of F_1 = F_2^gamma is replaced by the
comparison of the forward models: it fails after the convolution and at
low S/N (F119: 0.48 for a SiIV at S/N 50, before the fix of the frame of
7c.1), and it gives no credibility. It survives as the pre-filter of
section 3.1, with a low threshold (decided: Guido, 18 November 2026).

## 6. Ly-a only when asked; the assumed Ly-a and the interlopers

F51 and F93 item 3: today every region of the forest without a doublet is
labelled Ly-a with score 0, even when Ly-a is not asked. Proposed:
- the Ly-a is a hypothesis like the others, present only when asked
  (`Ly_a` or `Ly_ab` among the species); otherwise the region is
  *unidentified*;
- a region labelled Ly-a by default (no Ly-b covered, no alternative
  tested) is marked *assumed*, one with its Ly-b tested *identified*; the
  output, the window of the identifications and the log say which
  (F93: "the assumed Ly_a as such");
- the interlopers: in the forest, a metal candidate is compared with the
  Ly-a lines that would explain the same regions (section 3.2), with the
  Ly-b of the Ly-a hypothesis where covered (`doublet_vs_lyab` of the
  experiment); on the red side the alternative of a metal is the
  unidentified line;
- the forest treated statistically (the first point of the conversation
  of PR #798): if 7b adopts the mean transmission expected at each
  redshift as the backbone of the continuum, the Ly-a hypothesis of a
  group of the forest can become a statistical one, the likelihood of
  the flux under the distribution of the forest at that z, instead of a
  list of Voigt lines, with no free parameters per line. The metals, the
  strong H I and the DLAs stay explicit models. It removes the penalty of
  the free lines of the interlopers (section 4) and gives the prior of
  section 3.4 a basis in the data; an option to evaluate after 7b.1, not
  in the prototype. Likewise a learned proposer of hypotheses (the second
  point of #798, after phase 7) would only add hypotheses to a group,
  judged by the same likelihood and logged as the others.

## 7. Tellurics as a hypothesis

The identification reads the column `telluric` (the transmission, 1 where
there is none), whoever wrote it: 7b.6, or a model made outside Astrocook
(Molecfit, TelFit, the pipeline) imported as a column (after Guido's
reading, 18 November 2026). Three cases:
- *a model in the column, the flux not corrected*: as below, every
  hypothesis sees cont × telluric × absorption;
- *the flux already corrected* (the tellurics divided out, outside or by
  a recipe): the spectrum says so (a flag of its metadata, e.g.
  `TELLURIC_CORRECTED`, set by the recipe or the user), and the
  identification adds no telluric hypothesis; where the model is deep
  (transmission below about 0.5, when the column is there) the correction
  leaves residuals and inflated errors, and those regions are labelled
  *telluric residual*, kept out of the anchors (a doublet on a residual of
  the A band is the F54 case again);
- *nothing known*: the mask of the bands from the template of 7a.1, below.

F54 and F155: the O2 A band gives doublets with good scores (CIV at
z = 3.927 on J0103-1305), and the robustness grid has 172 false candidates
on tellurics (MgII 57 of its 64 false ones). Proposed:
- with the telluric model of 7b.6, the spectrum seen by every hypothesis
  is cont × telluric × absorption: all hypotheses include the telluric
  transmission (a column of the spectrum, last section), and the
  *telluric* hypothesis is the one with no absorber in the region, k = 0
  locally (the parameters of the telluric fit are global, a few for the
  whole spectrum); a doublet on a telluric line then has to beat the
  telluric alone with its 3 extra parameters;
- before 7b.6, or without the model: a mask of the bands from the template
  of 7a.1 (transmission below 0.9 at airmass 1.5), whose regions are
  labelled *telluric (masked)* and given no metal candidate; a weaker
  rule, measured against the model on the `telluric` case of the
  robustness grid;
- the application is 7c.5 (after 7b.6).

## 8. The standard multiplets (#455)

Today `STANDARD_MULTIPLETS` is the list of the lines both searched and
fitted, and only the doublets are searched. Proposed: two roles apart.
- *Search*: the anchor pairs, the two lines of a species most likely to be
  detected together (the strongest, or the doublet), as today for CIV,
  SiIV, MgII, plus pairs for FeII (2382/2600, 2344/2586) and SiII
  (1260/1526, 1304/1526) when they are asked; an anchor gives a kinematic
  candidate.
- *Model*: the hypothesis of the species at z has all its lines in the
  coverage of the spectrum, from the table of transitions, those weaker
  than a threshold at the column of the candidate dropped (their
  absorption below, say, 0.1 sigma); so the other lines are evidence
  (or counter-evidence: a FeII 2382 without its 2600 is not FeII).
- The 10 species of #455 with 3 or more transitions and no entry (in
  the table: CI, CoII, FeI, MgI, OI, PII, SIII, TiII, and NII, OII, OIII
  with their emission lines) get their entry of the model role
  automatically from the table; no new anchors for them until an
  analysis needs them.
- The fit of a subset of the lines of a multiplet (the request of #455)
  is the fit side, not the identification: the model role gives the
  lines, the fit window the subset (7d.3, F148).

## 9. Pairs and alternatives kept: requirements of 7c.3

The data structure of 7c.3 must hold, per group:
- the regions (x ranges, the detection significance) and the edges that
  made the group (the transition pair, the lag, the pre-filter value):
  the pairs of F93 item 1, stored, not found again by the window;
- per hypothesis: the components (series, z, logN, b, their errors, the
  uuids of the system list when it is the accepted one), the stage (AOD
  or fit), chi2, k, n, every criterion (AIC, AICc, BIC, SpIC H/a/b), the
  weight, the status (accepted, alternative, dropped by the screen), the
  reason in words for the log;
- the dropped candidates with their score and threshold (F119, today in
  `details.dropped_candidates`);
- the parameters of the run (criterion, prior, thresholds) to replay it.
It must live in the archive next to the system list (FITS cannot hold it
in a header: the mechanism shared with the Bayesian samples, F28), and be
read by the GUI (7c.7) and the MCP outputs (7c.6). Size: a few hundred
groups of a few hypotheses each, kB to a few hundred kB per spectrum.

Implemented in PR 7c.3: `astrocook/core/identification.py` (the structure,
in the metadata of the spectrum) and `astrocook/io/archive_data.py` (its
file in the archive, with the Bayesian samples).

## 10. The cost

**Measured.** `evaluate_model` (PR 7a.5): 0.7 s on a whole mock of 58
thousand pixels and 140 components (the model 0.05 s, the rest the
profiles one by one), 0.07 s on a 5 nm window. In the experiment (windows
of a few hundred pixels, one to two components per hypothesis): a median of 12-40 ms per trial for two hypotheses, each evaluated by AOD,
fitted and evaluated again, on mocks of a few thousand pixels; on the mocks
of `doublet_vs_lyab`, which span 420-505 nm (some 30 thousand pixels at R
45000) for 360 pixels evaluated, 1.05 s (p95 4 s): `evaluate_model` and the
fit compute the model and the profiles on the whole spectrum, so the cost
follows the length of the spectrum, not the pixels of the comparison.

**Estimated per spectrum.** A mock of the `full` grid at z_em 3.5 has
about 1600 forest lines (12 < logN < 17) in its forest, a few hundred
regions after the merging, and some tens on the red side. With n + 1
hypotheses per group, about 1000-2000 hypotheses per spectrum; at 10-50 ms
each (AOD screen and the fit of the finalists of groups of a few
components) that is 0.5-2 minutes, against the hour of section 5. Two
conditions for it:
- the evaluation in a window around the lines of the group, not on the
  whole spectrum: the window of the pixels of the group (±200 km/s around
  each line in the experiment; the union of ±(3 FWHM + 5 b)
  around each line, decided), and the profiles of
  `evaluate_model` computed on that slice;
- the components outside the group fixed: their transmission computed
  once per group and multiplied (a cache, not a recomputation per
  hypothesis); the system list of the rest does not change between the
  alternatives of a group.
The groups are independent, so the evaluation can be parallel; the budget
is a time per spectrum logged next to the result.

**Done in PR 7c.4a** (18 November 2026): `EvaluationWindow`
(`astrocook/fitting/evaluation_window.py`) cuts the pixels of a group out of
the spectrum and keeps the fixed components that reach them, once per group;
the fixed components are frozen in the model of the window rather than
multiplied as a transmission (the product of the convolved fluxes moved the
ΔIC of blended hypotheses by up to 1.3 at R 8900); the short fit takes its
pixels from the data as on the whole spectrum, not from the window of the
start (limited to it, a broad saturated line started narrow was fitted on
its core); the pixels of the comparison are the union of the windows and
of the pixels each hypothesis takes from the data (Guido, 18 November 2026).
The statistics equal those on the whole spectrum within 1e-5; a hypothesis takes 17-22 ms in a
forest of 32 thousand pixels, the trials of `doublet_vs_lyab` 1.34 s → 0.068 s
(the row of 7c.4a in `ROADMAP.md`).

## 11. The criterion of "better" for 7c.8

Measured with `validate_mocks.py --compare` on both grids, the pairs with
`validate_lyab_pairs.py`, and the development spectra (7a.4).

Must improve:
- the completeness of the metals between 3 and 10 sigma (F57: CIV 0.12 and
  0.55, SiIV 0.17 and 0.45, MgII 0.29 at 5-10), in both grids;
- the pixel outcomes of 7a.2: the detectable absorbed pixels with a tested
  identification (about 1% today) or an explicit unidentified label; the
  pixels in no region (10-22%);
- the false candidates on tellurics (172 in `robustness`) and on the
  associated absorbers (49 CIV, 39 SiIV);
- the completeness of the Ly-a/Ly-b pairs (0.087), at a purity of at
  least 0.84.

Must not worsen (beyond the noise of the seeds):
- the purity of CIV, SiIV, MgII (`full` 0.987-0.993, `robustness`
  0.82-0.86);
- the fit after the identification: matched, spurious, the pull widths
  and the system column of `doublets_auto`;
- the time: within the budget of section 10, and stated.

New, measured by 7c.4 (F164, proposed in PR #798): the calibration of
the credibility. Of the hypotheses given weight 0.8, about 80% should be
right: a reliability diagram (bins of weight against the fraction of
hypotheses that are true on the mocks), its expected calibration error
and the Brier score, in `validate_mocks.py`, over all the hypotheses of
the groups (the alternatives too, not only the accepted ones), per
category (metal against interloper, coincidence, telluric) and by S/N
and R; the same metric then for the flags of the components (7d). A
criterion of "better" in its own right: a prototype more complete but
overconfident is not better. On the mocks the noise is independent; on
the twins and the development spectra the calibration measures also the
tempering of section 3.3.

## 12. The interface of the prototype

- `absorbers.identify_lines` gets `method`, `"doublets"` (the current
  algorithm, the default) or `"hypotheses"` (the prototype), as the
  `method` of the continuum of 7b.2; the parameters of the prototype:
  `species` (the anchors searched; default all the species with an anchor,
  below), `criterion` (`"spic"`, `"bic"`,
  `"aic"`), `prior` (`"none"`, `"incidence"`), `assume_lya` (the Ly-a by
  default of the forest, off unless asked, F51), `max_group` and the screen
  threshold; the rest of the chain (`doublets_auto`, `components_auto`)
  reads the same `abs_ids` and system list, so it works unchanged with
  either method.
- *The species by default* (after Guido's question, 18 November 2026: why
  ask the user at all?). Proposed: all, in the sense of section 8 and 3.5:
  every species with an anchor is searched (CIV, SiIV, MgII, NV, OVI,
  AlIII, CII, FeII, SiII, ... and `Ly_ab` in the forest), every
  species of the table is swept at the z of the anchored systems; the
  parameter restricts the anchors, it does not need to be set. The cost of
  "all" is not the time but the purity: each anchor searched adds chances
  of a coincidence (the hard case of section 4), the look-elsewhere effect
  of section 5 grows with the number of species; the prior of the
  incidence (section 3.4, an option) is the natural guard for the rare
  species. 7c.4 measures the purity with all the anchors against CIV,
  SiIV, MgII alone; if it drops, the default is the common anchors, the
  rest swept.
- *The fit of the identification* (after Guido's question). The short fit
  of the finalists consolidates the comparison; it is not the model. The
  accepted hypothesis is the output system list, with the parameters of
  its short fit as the starting values (better than the AOD guesses of
  today) and marked as such; the alternatives keep theirs in the
  structure of 7c.3, so that the GUI can show them and the user swap one
  in. The model is made by the fits that follow (`doublets_auto`,
  `components_auto`: the components added from the residuals, the full
  number of iterations, the flags of 7d).
- *Showing the graph* (after Guido's question; the GUI is 7c.7). The
  graph as such (nodes and edges) says little; on the spectrum it does:
  the regions coloured by their accepted label, their transparency by its
  weight, and the edges of the selected region drawn as arcs to the
  regions its hypotheses link (one colour per species), so that a click on
  a line shows which lines its identification rests on. For a group, its
  hypotheses side by side as velocity stacks of their lines (the panel of
  the System Inspector and of `spectrum_plot_velocity`), with their
  weights. Those panels should come from one piece of code: F163 (the
  Inspector draws its panel apart from `plotting.py`) is best closed
  before or with 7c.7, so that the GUI and the MCP tools show the
  alternatives in the same way.
- *Replacing the old method* (after Guido's question). The prototype
  becomes the default only in 7c.8, if the benchmark shows it better, and
  the old method stays selectable (`method="doublets"`) at least for a
  release; whether to remove it then is Guido's decision. One care for
  7c.8: the logs record the parameters used, and a log written before the
  change has no `method`; it must replay with `"doublets"`, not with the
  new default (a default by version of the log, or `method` added on
  load). the system list of the accepted hypotheses as today, the
  groups with their alternatives in the structure of 7c.3, a summary in
  `RecipeResult.extra` (7c.6), and one line of the log per group decided
  ("CIV at z = 2.2295 (w = 0.97) against two Ly-a (0.03), SpIC").
- The system of an accepted Ly-a/Ly-b pair (from 7c.1): an H I system
  with its Ly-b, one component fitted on both lines, and the Ly-b region
  labelled as such (decided `Ly_ab`: Guido, 18 November 2026); a `Ly_a`
  system where the Ly-b is not covered. The higher orders of the series
  (Ly-g and beyond) are not used by the identification, neither as
  anchors nor in the model (Guido, 18 November 2026, after a first
  proposal of the whole `Lyman` series):
  - they fall in ever denser absorption (the Ly-b, Ly-g, ... forests
    overlap) and converge, so in the identification they would mostly
    find absorption to explain whatever the hypothesis: the pixels of the
    comparison hold unmodelled blends, and a hypothesis that predicts a
    line there is rewarded for a blend it did not cause. A line predicted
    on a clean continuum can falsify an H I hypothesis, but a line found
    in a blend supports nothing; the Ly-b already has this bias in the
    Ly-b forest (the `doublet_vs_lyab` of section 4 is its best case);
  - they are not always accessible: out of the spectrum, or below the
    Lyman limit of a stronger system at higher z, where the flux is zero;
  - they are hard to fit even when identified, and that belongs to the
    fit: the column of a saturated H I (logN above about 15, where Ly-a
    and Ly-b are flat-bottomed) needs them, as a choice of the lines to
    fit (#455, 7d), not of the identification.
  The software can carry the whole series (checked: a `Lyman` and a
  `Ly_ab` system on a mock that covers Ly-a and Ly-b only give the same
  logN, b and chi2, the lines out of the spectrum add nothing), so the
  choice can be reopened for the fit of the strong systems. The
  Lyman-limit continuum belongs to the DLAs and LLSs of 7b.

## Common decisions with 7b

Taken to the conversation "7 · plan of the phase" (Guido, 16 November
2026); the three proposals below accepted by Guido on 18 November 2026:
1. *One vocabulary of the regions*: BAL, DLA, telluric, forest (Ly-a,
   assumed or identified), metal, unidentified, as the categories of the
   truth map of 7a.1 (forest, metal, strong_hi, bal, telluric). 7b writes
   the regions it owns (the masks of the BAL troughs and of the DLA cores
   and wings, the `telluric` column), the identification reads them, takes
   them out of its groups, and writes the final label of every region in
   its structure (7c.3) and in `abs_ids`.
2. *The continuum of 7c.4*: the current one (`estimate_auto`), so that the
   two prototypes are measured apart; 7c.4 measured again on the
   continuum of 7b.2 when it exists.
3. *The telluric model as a column*: `telluric`, the model transmission
   (1 where there is none), written by 7b.6 or imported (section 7), read
   by the identification and the fit.

## Decisions (Guido, 18 November 2026)

On the questions of the first version of this document:
1. The criterion: SpIC_H by default, BIC reported next to it.
2. The prior: none in the prototype (an option of the recipe, with its
   sources), the metals of the forest off by default and opened by a link
   to another ion.
3. The detection by matched filters, threshold 3 sigma per line and on
   the joint significance of a doublet, as the direction of 7c.4.
4. The R^2 score kept as the pre-filter of the candidates.
5. The groups of the identification, distinct from those of the fits
   (section 3.1): at most 6 regions, n + 1 hypotheses.
6. The window of the pixels of a group: ±(3 FWHM + 5 b) around each line.
7. An accepted Ly-a/Ly-b pair gives an H I system with its Ly-b (`Ly_ab`;
   `Ly_a` where the Ly-b is not covered); the higher orders not used by
   the identification (section 12).
8. The two roles of the multiplets (search anchors and model lines), with
   the anchors of FeII and SiII.
9. The three common decisions with 7b, above.

Added after the reading, and accepted by Guido the same day: the sweep of the anchored
systems for the rarer and single lines (section 3.5); the species "all"
by default (section 12); the short fit as starting values, not the model
(section 12); the display of the graph and F163 (section 12); the
telluric model made outside Astrocook, and the flux already corrected
(section 7); the replay of the old logs after 7c.8 (section 12); the
calibration of the credibility, F164, as a criterion of "better"
(section 11), with the weights tempered for a chi2 that is not calibrated
(section 3.3); the forest as a statistical hypothesis, an option after
7b.1 (section 6).
