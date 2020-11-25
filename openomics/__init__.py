from __future__ import print_function, division, absolute_import
import pandas as pd
import dask.dataframe as dd

# -*- coding: utf-8 -*-

"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = 'nhat.tran@mavs.uta.edu'
__version__ = '0.7.9'

__BACKEND__ = "pandas"
backend = pd

from . import database

from .transcriptomics import (
    ExpressionData, MessengerRNA, MicroRNA, LncRNA,
)

from .genomics import (
    SomaticMutation, DNAMethylation, CopyNumberVariation
)

from .proteomics import (
    Protein
)

from .clinical import ClinicalData

from .multiomics import (
    MultiOmics
)

from .visualization import (
    umap
)


def set_backend(new):
    global __BACKEND__
    global backend
    assert new in ["dask", "pandas"]

    __BACKEND__ = new
    if new == "dask":
        backend = dd
    else:
        backend = pd
