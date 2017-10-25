# Astrocook

A thousand ways to cook a spectrum!

## Getting Started

To get a copy of Astrocook on your local machine: 

```
git clone https://github.com/DAS-OATs/astrocook
```

### Prerequisites

Astrocook requires the following packages to run:

* [Astropy](http://www.astropy.org), including [Specutils](http://specutils.readthedocs.io/en/latest/)
* [SciPy](https://www.scipy.org), and in particular [NumPy](http://www.numpy.org) and [matplotlib](https://matplotlib.org) 
* [LmFit](https://lmfit.github.io/lmfit-py/)
* [StatsModels](http://www.statsmodels.org/stable/index.html)
* [Cycler](https://pypi.python.org/pypi/Cycler)


## Running the tests

The following tests are available:

* test.py: import and plot a spectrum, create a power-law model, and deredden the spectrum.
* line_test.py: create a list of absoprtion lines from a spectrum and fit them with Voigt profiles.
* syst_test.py: create a list of CIV doublets from a spectrum and fit them with Voigt profiles.

To run the tests:

```
python <name-of-the-test>
```

DISCLAIMER: Due to the ongoing development, some tests may occasionally crash.

## Contributing

A CONTRIBUTING.md file will be soon uploaded to detail our code of conduct and the process for submitting pull requests to us.

## Authors

* **[Guido Cupani](https://github.com/gcupani)** - [INAF-OATs](http://www.oats.inaf.it/index.php/en/)
* **[Giorgio Calderone](https://github.com/gcalderone)** - [INAF-OATs](http://www.oats.inaf.it/index.php/en/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Acknowledgments

The project is carried out at [INAF-OATs](http://www.oats.inaf.it/index.php/en/) with the contributions of Stefano Cristiani and Giuliano Taffoni.