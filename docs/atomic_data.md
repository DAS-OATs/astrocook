# Atomic Data

Astrocook takes the rest wavelengths, oscillator strengths, damping
constants and masses of the transitions from a single table,
`astrocook/data/atomic_data.ecsv`. The fitter, the identification, the
mock spectra and the GUI (redshift cursor, line markers, System Inspector)
all use this table. This page describes its content and where each value
comes from.

## The table

The table is an [ECSV](https://docs.astropy.org/en/stable/io/ascii/ecsv.html)
file, so astropy can read it with its units and the list of references:

```python
from astropy.table import Table
from astrocook.core.atomic_data import ATOMIC_DATA_PATH

table = Table.read(ATOMIC_DATA_PATH, format="ascii.ecsv")
print(table.meta["references"])
```

Inside Astrocook the same values are available as `ATOM_DATA` (a dictionary
keyed by transition) and `xem_d` (rest wavelengths as astropy quantities),
both in `astrocook.core.atomic_data`.

| Column | Unit | Content |
|---|---|---|
| `name` | | transition: ion and rounded wavelength, e.g. `CIV_1548`, `Ly_a` |
| `wave` | Å | rest wavelength in vacuum |
| `f` | | oscillator strength (0 for emission-only lines and edges) |
| `gamma` | s⁻¹ | damping constant |
| `mass` | u | atomic mass |
| `wave_ref`, `f_ref`, `gamma_ref` | | source of `wave`, `f` and `gamma` (see below) |

The table has 211 transitions: the H I Lyman series up to Ly-15 and the
Lyman limit, D I, the metal absorption lines of the intergalactic and
interstellar medium, and some emission lines (used as markers, with
`f = 0`).

A source key means that the value was taken from that source, or agrees
with it within 2 mÅ in wavelength, 1% in `f` and 2% in `gamma`.

| Key | Source | Used for |
|---|---|---|
| `M03` | Morton D. C. 2003, ApJS 149, 205 ([ADS](https://ui.adsabs.harvard.edu/abs/2003ApJS..149..205M)) | most of the table: 159 wavelengths, 150 `f`, 129 `gamma` |
| `K14` | Kisielius R. et al. 2014, ApJ 780, 76 ([ADS](https://ui.adsabs.harvard.edu/abs/2014ApJ...780...76K)) | `f` of S II 1250, 1253, 1259 |
| `K15` | Kisielius R. et al. 2015, ApJ 804, 76 ([ADS](https://ui.adsabs.harvard.edu/abs/2015ApJ...804...76K)) | `f` of Zn II 2026, 2062 |
| `BB19` | Boissé P., Bergeron J. 2019, A&A 622, A140 ([ADS](https://ui.adsabs.harvard.edu/abs/2019A&A...622A.140B)) | `f` of Ni II 1317, 1370, 1454, 1709, 1741, 1751 |
| `PA04` | Petitjean P., Aracil B. 2004, A&A 422, 523 ([ADS](https://ui.adsabs.harvard.edu/abs/2004A&A...422..523P)) | wavelengths of C IV 1548, 1550 |
| `VPFIT` | the `atom.dat` file of VPFIT (copy received from P. Noterdaeme, October 2026), which cites its sources line by line | values that Morton (2003) does not list: C II* 1335, Si II* 1264, O II 832-834, O III 832, Ca I 4227, Ca II 3934, 3969, Fe I 1851, Ni II 1317 (wavelength), and 29 `gamma` |
| `V1` | `atom_par.dat` of Astrocook 1; original source not recorded | Balmer lines, emission lines, some `gamma`, C II 903, 904, Fe III 859, Zn I 1589, Sn II 1400 |
| `V2` | the table of Astrocook 2.0-2.1; original source not recorded | the Lyman limit, [Ne III] 3869 and [Ne V] 3346 (vacuum wavelengths) |
| `-` | not applicable | values equal to 0 |

The masses are standard atomic weights, carried over from the earlier
tables without a recorded source.

## Choices and differences from other line lists

The table was rebuilt in 2026 (issue
[#568](https://github.com/DAS-OATs/astrocook/issues/568)). Each value was
compared with Morton (2003), as distributed with
[linetools](https://github.com/linetools/linetools), and with the `atom.dat`
of VPFIT. Some values were also compared with the line list of
[VoigtFit](https://github.com/jkrogager/VoigtFit).

**Newer values adopted instead of Morton (2003).** They follow the `atom.dat`
of VPFIT, and each one is checked against its original paper. The column
"Before" gives the values of Astrocook 2.1, which are those of Morton
(2003) except Ni II 1317, a line that Morton does not list:

| Lines | Before | Adopted | Effect on the column density |
|---|---|---|---|
| Zn II 2026, 2062 | f = 0.501, 0.246 | f = 0.630, 0.309 (K15) | logN(Zn II) 0.10 dex lower |
| Ni II 1317, 1370, 1454 | f = 0.0779, 0.0769, 0.0323 | f = 0.0596, 0.0616, 0.0220 (BB19) | logN(Ni II) up to 0.17 dex higher |
| Ni II 1709, 1741, 1751 | f = 0.0324, 0.0427, 0.0277 | f = 0.0318, 0.0432, 0.0280 (BB19) | within 0.01 dex |
| S II 1250, 1253, 1259 | f = 0.00543, 0.0109, 0.0166 | f = 0.00515, 0.01035, 0.015675 (K14) | logN(S II) 0.02 dex higher |
| C IV 1548, 1550 | 1548.204, 1550.781 Å | 1548.2049, 1550.7785 Å (PA04) | none; the lines move by 0.2 and 0.5 km/s |

Petitjean & Aracil (2004) measured the C IV wavelengths from high-redshift
quasar absorbers, calibrated on Si IV. Their values agree with the
laboratory wavelengths of Griesmann & Kling (2000, ApJ 536, L113,
[ADS](https://ui.adsabs.harvard.edu/abs/2000ApJ...536L.113G)).

**Values of VPFIT not adopted**, because their source is not clear:

- the `gamma` of O I 971, 976, Mg I 2026, Mg II 2796, 2803, Co II 1466 and
  Fe I 3441, which differ from those of the table;
- the placeholder `gamma` of 2e8 s⁻¹ for O II, O III and Fe II 1081;
- S II 1253 at 1253.811 Å (Morton: 1253.805 Å);
- the `f` of N I 1134.4149 (0.0278; Morton: 0.0287);
- the values of Sn II 1400.

These choices only matter for the damping wings of very high columns, or
for a wavelength shift of 1.4 km/s on S II 1253.

**Corrections to earlier versions** (Astrocook 2.1 and earlier):

- The wavelengths of O III 832, P III 913, P V 1117 and Fe I 2501 were
  wrong by 0.2-0.7 Å, from a slip in the conversion from nm to Å.
- The damping constant of Ar I 1066 was ten times too small.
- The wavelength of Al II 1670 was 0.2 km/s off.
- The damping constant of C I 1560 was 1.17e8 s⁻¹ instead of 1.27e8.
- [Ne III] 3869 and [Ne V] 3346 were given in air instead of vacuum.
- The table listed a Ti II line at 3066.35 Å that no reference lists. It is
  kept with `f = 0`, so that sessions naming it still open, and the real
  Ti II 3067.24 Å line is added.

**Transitions of Astrocook 1 not in the table.** The V1 file also listed H₂
and CO bands, the Paschen and Brackett series, Si I, Balmer lines beyond
Hε, and stellar and blended emission markers. Astrocook 2 never used them.
They can be added when an analysis needs them.

## Series and multiplets

A series is written as a list of transitions (`CIV_1548,CIV_1550`) or by a
name. Astrocook knows two kinds of names:

- **Standard multiplets** (`STANDARD_MULTIPLETS`): the groups used by the
  identification and the automatic fits, strongest line first.
  - Doublets: `CIV`, `SiIV`, `NV`, `OVI`, `MgII`, `AlIII`, `CII` (with
    C II* 1335), `ZnII`, `NaI`, `CaII`.
  - Multiplets: `NI`, `SiII`, `SII`, `CrII`, `MnII`, `FeII`, `NiII`.
  - Hydrogen: `Lyman` (Ly-α to Ly-8), `Ly_a`, `Ly_ab`.
- **Series of the V1 syntax**, used by the redshift cursor of the GUI and
  for older sessions. Each is the set of the transitions in the table with
  the same name before the underscore (e.g. `MgII` gives its six lines from
  1025 to 2803 Å, `Ly` the whole Lyman series), in the order of the table.
  Their first transition is the reference of the redshift cursor: Ly-α for
  `Ly`, `D` and `H`. The aliases `Ly-a`, `Ly-ab` and `Ly-abg` give the first
  one, two or three Lyman lines.

A whole name is looked up first among the standard multiplets: `MgII` alone
is the doublet, while in a list such as `MgII;Ly_b` it is the V1 series.
