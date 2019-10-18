# Astrocook

A thousand ways to cook a spectrum!

 [![DOI](https://zenodo.org/badge/78840469.svg)](https://zenodo.org/badge/latestdoi/78840469)

## Getting Started

To get a copy of Astrocook on your local machine:

```
git clone https://github.com/DAS-OATs/astrocook
```

### Prerequisites

Astrocook requires the following packages to run:

* [Astropy 3.2.2](http://www.astropy.org)
* [SciPy 1.3.1](https://www.scipy.org)
* [Specutils 0.6](http://specutils.readthedocs.io/en/latest/)
* [NumPy 1.17.3](http://www.numpy.org)
* [LmFit 0.9.14](https://lmfit.github.io/lmfit-py/)
* [Cycler 0.10.0](https://pypi.python.org/pypi/Cycler)
* [StatsModels 0.10.1](http://www.statsmodels.org/stable/index.html)
* [matplotlib 3.1.1](https://matplotlib.org)
* [Sphinx 2.2.0](http://www.sphinx-doc.org/en/master/)


## Running the software

The Astrocook user manual is under construction. In the meanwhile, you can try the following commands:

* ```python mock_demo.py```, which creates a mock spectrum, plots it and saves it into an Astrocook archive (.acs).
* ```python ac_gui.py```, which launches the (self-explaining) Astrocook GUI to perform the analysis.

A note on .acs archive: they are normal tarballs containing all the products of an analysis session (spectrum, optional list of lines and absorption systems, etc.) in FITS format. To extract an .acs archive: ```tar -zxvf [name].acs```.

**More information coming soon!**


## Contributing

A CONTRIBUTING.md file will be soon uploaded to detail our code of conduct and the process for submitting pull requests to us.

## Authors

* **[Guido Cupani](https://github.com/gcupani)** - [INAF-OATs](http://www.oats.inaf.it/index.php/en/)
* **[Giorgio Calderone](https://github.com/gcalderone)** - [INAF-OATs](http://www.oats.inaf.it/index.php/en/)
* **[Stefano Alberto Russo](https://github.com/sarusso)** - [INAF-OATs](http://www.oats.inaf.it/index.php/en/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Releases

v0.1 - Start of the project
v0.2 - First release with a DOI; GUI available; tools for finding lines, determining continuum, and fitting systems
v0.3 – Transition to Python 3; complete makeover of classes and methods

## License

Astrocook is licensed under the [GNU General Public License (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Acknowledgments

The project is carried out at [INAF-OATs](http://www.oats.inaf.it/index.php/en/) with the contributions of Stefano Cristiani and Giuliano Taffoni.
