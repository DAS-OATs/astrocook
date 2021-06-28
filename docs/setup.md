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


## Install the software

You can download the latest release of Astrocook (**v1.0.0-rc.1**) from [this page](https://github.com/DAS-OATs/astrocook/releases/tag/v1.0.0-rc.1). Alternatively, you can clone the GitHub repository on your local machine and checkout the specific tag:
```
$ cd /your/path/
$ git clone https://github.com/DAS-OATs/astrocook
$ git checkout tag/v1.0.0-rc.1
```

Some additional features may be available in the `develop` branch before they are included in an official release. To fetch the last commit of the `develop` branch:
```
$ cd /your/path/astrocook/
$ git fetch origin develop
$ git checkout develop
```

## Troubleshooting

Occasionally, Astrocook may behave erratically. If this happens, you are strongly encouraged to [report the bug](mailto:guido.cupani@inaf.it).

❗️ **If Astrocook stops responding, you can kill the GUI with `ctrl+C`. Since this will destroy all sessions, you are suggested to frequently save your analysis.**
