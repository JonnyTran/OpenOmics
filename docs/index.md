# Welcome to OpenOmics's documentation!

[![PyPI version](https://badge.fury.io/py/openomics.svg)](https://badge.fury.io/py/openomics)
[![Documentation Status](https://readthedocs.org/projects/openomics/badge/?version=latest)](https://openomics.readthedocs.io/en/latest/?badge=latest)
[![OpenOmics](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics/branch/master/graph/badge.svg?token=WAN3PJwM42)](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics)

This Python package provide a series of tools to integrate and query the genomics, transcriptomics, proteomics, and
clinical data (aka, the multi-omics data). With scalable data-frame manipulation tools, OpenOmics facilitates the common
coding tasks when preparing data for RNA-seq bioinformatics analysis.

## Features

- Provides a bioinformatics workflow to generate integrative results from multi-omics data.
- Facilitates integration of various bio-databases, multi-omics expression, genomics, and clinical data.
- Highly flexible to different data types and missing data.
- Provides researchers with means to consistently store and explore their experimental datasets.
- Enables scalable performance with parallel computing, while easily configurable to deploy on both single machine and a
  cluster.

## Table of Content

```{toctree}
:maxdepth: 2
:caption: Using OpenOmics Python API
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
:caption: API Reference

modules/openomics.multiomics
modules/openomics.annotate
modules/openomics.database.annotation
modules/openomics.database.sequence
modules/openomics.database.interaction
modules/openomics.database.disease
modules/openomics.database.ontology
modules/openomics.utils
```

```{toctree}
:maxdepth: 1
:caption: MISC

misc/faq
```

```{toctree}
:maxdepth: 1
:caption: Contributing and releases
contributing
history
```
