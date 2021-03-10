# Welcome to OpenOmics's documentation!

![Logo](https://github.com/BioMeCIS-Lab/OpenOmics/raw/master/openomics_web/assets/openomics_logo.png)

[![PyPI version](https://badge.fury.io/py/openomics.svg)](https://badge.fury.io/py/openomics)
[![Documentation Status](https://readthedocs.org/projects/openomics/badge/?version=latest)](https://openomics.readthedocs.io/en/latest/?badge=latest)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![OpenOmics](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/JonnyTran/OpenOmics/branch/master/graph/badge.svg)](https://codecov.io/gh/JonnyTran/OpenOmics)

This Python package provide a series of tools to integrate and query the genomics, transcriptomics, proteomics, and
clinical data (aka, the multi-omics data). With scalable data-frame manipulation tools, OpenOmics facilitates the common
coding tasks when preparing data for RNA-seq bioinformatics analysis.

- Free software: MIT license
- Documentation: https://openomics.readthedocs.io.

## Features

- Provides a bioinformatics workflow to generate integrative results from multi-omics data.
- Facilitates integration of various bio-databases, multi-omics expression, genomics, and clinical data.
- Highly flexible to different data types and missing data.
- Provides researchers with means to consistently store and explore their experimental datasets.
- Enables scalable performance with parallel computing, while easily configurable to deploy on both single machine and a
  cluster.

## Table of Content

```{toctree}
:maxdepth: 1
:caption: Using OpenOmics
:name: mastertoc
installation
usage/getting-started
usage/import-your-dataset
usage/annotate-external-databases
usage/network-data
usage/preprocess-downstream-analysis
```

```{toctree}
:maxdepth: 1
:caption: Reference and Contributing
contributing
```

## Modules

```{autosummary}
:toctree: _autosummary
:template: custom-module-template.rst

openomics.multiomics
openomics.clinical
openomics.transcriptomics
openomics.database
```

```{eval-rst}
.. automodapi::

openomics
```

```{eval-rst}
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```
