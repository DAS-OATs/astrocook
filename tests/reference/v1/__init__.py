"""V1 reference for the V1/V2 parity tests (phase 3, ``ROADMAP.md``).

The V1 code is the validated reference behind the published papers. These
modules are verbatim copies of ``astrocook/legacy/`` at commit c5c7a86a
(October 2026), limited to what the parity tests use and its imports:
``spectrum`` (``Spectrum``), ``functions`` (``lines_voigt``,
``convolve_simple``, ``trans_parse``, ``x_convert``), ``vars``,
``cookbook_general``, ``cookbook_continuum`` (with
``cookbook_continuum_old``), ``frame``, ``line_list``, ``format`` and
``syst_list`` (used by ``Spectrum._rebin``), ``message`` and ``utils``. Import them as ``reference.v1.<module>``.

Only ``vars.py`` is changed (marked in the file): it reads the V1 data
files from ``data/`` (copies of ``atom_par.dat``, ``qso_composite.dat``
and ``table2_inoue.dat``), so that the reference does not change when the
data of the package do, and it does not load ``sky_telluric.fits``.

Do not edit these modules, except to keep them importable: a parity test
that fails must be explained by V2, not fixed here. A few code paths
import V1 modules that are not copied (e.g. the V1 session); the tests
do not reach them. Not linted (``pyproject.toml``).
"""
