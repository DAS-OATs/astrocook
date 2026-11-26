<img width="400" height="311" alt="logo_3d_HR" src="https://github.com/user-attachments/assets/89f51f8b-128a-4d06-982a-e925d24c43da" />

# Astrocook

*A thousand recipes to cook a spectrum.*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1346736.svg)](https://doi.org/10.5281/zenodo.1346736)

Astrocook is a Python package and graphical interface for the analysis of quasar spectra: estimate of the continuum, identification of absorption systems and fit of Voigt profiles, interactively and in a reproducible way. Every operation is a *recipe* with its own parameters, available from the menus, from Python scripts and from AI agents (through an MCP server); the history of the session is recorded and can be exported as a description of the methods.

## Installation

Astrocook needs Python 3.12 or later (the installation downloads PySide6, about 200 MB):
```
pip install astrocook
```
Optional extras: `astrocook[bayes]` (MCMC and nested sampling) and `astrocook[mcp]` (server for AI agents). Then run `astrocook`. On macOS a standalone app is available from the [releases page](https://github.com/DAS-OATs/astrocook/releases/latest).

## Documentation

User guide, tutorials and reference: <https://das-oats.github.io/astrocook/>. Bugs and questions: [GitHub issues](https://github.com/DAS-OATs/astrocook/issues).

## Credits

© 2017-2026 INAF – Istituto Nazionale di Astrofisica. Version: 2.9.1

Developed by **Guido Cupani** at INAF–Osservatorio Astronomico di Trieste.

**Many thanks** to (in alphabetical order): Giorgio Calderone, Stefano Cristiani, Simona Di Stefano, Valentina D'Odorico, Francesco Guarneri, Elena Marcuzzo, Stefano Alberto Russo, Andrea Trost.

The development of Astrocook was assisted by the AI models Claude (Anthropic) and Gemini (Google).

Logo and icon by **[Claudia Cupani](https://claudiacupani.com)**.

Partly funded by the PRIN PNRR project “Next generation computing & data technologies to probe the cosmic metal content”.

If you use Astrocook in your research, please cite: **[Cupani et al. 2026](https://scixplorer.org/abs/2026SPIE14155E..1XC/abstract) (Proc. SPIE 14155, 141551X)**. To cite the software itself, use its DOI, [10.5281/zenodo.1346736](https://doi.org/10.5281/zenodo.1346736) (it always points to the latest version; the page of each version has its own DOI). The `CITATION.cff` file in the repository holds the same information.

## License
Astrocook is free software, released under the [GNU General Public License v3.0 or later](https://github.com/DAS-OATs/astrocook/blob/master/LICENSE) (`GPL-3.0-or-later`).

## Disclaimer
THE SOFTWARE ASTROCOOK IS PROVIDED ​“AS IS”. THE DEVELOPERS OF ASTROCOOK MAKE NO OTHER WARRANTIES, EXPRESS OR IMPLIED, AND HEREBY DISCLAIM ALL IMPLIED WARRANTIES, INCLUDING ANY WARRANTY OF MERCHANTABILITY AND WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.
