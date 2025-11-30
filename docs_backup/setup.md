---
layout: default
title: Setup
parent: Getting started
nav_order: 1
---

# Setup
{: .no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

## Dependencies

Astrocook is written in Python 3 and requires the following packages:

* [Python 3.12+](http://www.python.org)
* [wxPython 4.2+](https://wxpython.org/)
* [Astropy 6.1+](http://www.astropy.org)
* [tqdm 4.66+](https://github.com/tqdm/tqdm)
* [matplotlib 3.9+ (base)](https://matplotlib.org)
* [SciPy 1.14+](https://www.scipy.org)
* [Lmfit 1.3+](https://lmfit.github.io/lmfit-py/)
* [SciencePlots 2.1+](https://github.com/garrettj403/SciencePlots)

* [Specutils 0.6](http://specutils.readthedocs.io/en/latest/)
* [NumPy 1.17.3](http://www.numpy.org)
* [Cycler 0.10.0](https://pypi.python.org/pypi/Cycler)
* [StatsModels 0.10.1](http://www.statsmodels.org/stable/index.html)

You are suggested to manage the dependencies using [Conda](https://docs.conda.io/projects/conda/en/latest/).
The quickest way to do it is to get the Miniconda installer (instructions [here](https://docs.conda.io/en/latest/miniconda.html)) and to open an environment dedicated to Astrocook:
```
$ conda create -n astrocook python=3.12
$ source activate astrocook
```

This is an example of how you can install the required packages from the Conda repository:  
```
$ conda install conda-forge::wxpython
$ conda install conda-forge::astropy
$ conda install conda-forge::tqdm
$ conda install conda-forge::matplotlib-base
$ conda install conda-forge::scipy
$ conda install conda-forge::lmfit
$ conda install conda-forge::scienceplots

$ conda install -c anaconda package-name=package-version
```
If this doesn't work, look for the package in the [Anaconda Cloud](https://anaconda.org/).


## Install the package

1. Clone the GitHub repository on your local machine:
```
$ cd /your/path/
$ git clone https://github.com/DAS-OATs/astrocook
```
2. Fetch the last commit of the `develop` branch
```
$ cd /your/path/astrocook/
$ git fetch origin develop
$ git checkout develop
```

❗️ **The current documentation applies to v1.0.0, which has not been released yet. You are suggested to work on the `develop` branch until it gets merged on the `master` branch.**
