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

* [Astropy 4.0](http://www.astropy.org)
* [SciPy 1.3.1](https://www.scipy.org)
* [Specutils 0.6](http://specutils.readthedocs.io/en/latest/)
* [NumPy 1.17.3](http://www.numpy.org)
* [LmFit 1.0.0](https://lmfit.github.io/lmfit-py/)
* [Cycler 0.10.0](https://pypi.python.org/pypi/Cycler)
* [StatsModels 0.10.1](http://www.statsmodels.org/stable/index.html)
* [matplotlib 3.1.1](https://matplotlib.org)
* [Sphinx 2.2.0](http://www.sphinx-doc.org/en/master/)
* [wxPython 4.0.4](https://wxpython.org/)
* [tqdm 4.36.1](https://github.com/tqdm/tqdm)

You are suggested to manage the dependencies using [Conda](https://docs.conda.io/projects/conda/en/latest/).
The quickest way to do it is to get the Miniconda installer (instructions [here](https://docs.conda.io/en/latest/miniconda.html)) and to open an environment dedicated to Astrocook:
```
$ conda create -n astrocook python=3.7
$ source activate astrocook
```

This is an example of how you can install the required packages from the Conda repository:  
```
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
