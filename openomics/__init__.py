from __future__ import print_function, division, absolute_import

import astropy
import logging
import os
import sys

import dask.dataframe as dd
import pandas as pd
"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = "nhat.tran@mavs.uta.edu"
__version__ = "0.8.5"

# Initialize configurations
this = sys.modules[__name__]
this.config = {}

# Read configuration from ~/.openomics/conf.json

# Set pandas backend
if not this.config:
    this.config["backend"] = pd

    # Set cache download directory
    this.config["cache_dir"] = astropy.config.get_cache_dir(this.__name__)
    logging.info("Cache directory is", this.config["cache_dir"])

from . import database, utils

from .transcriptomics import (
    Expression,
    MessengerRNA,
    MicroRNA,
    LncRNA,
)

from .genomics import SomaticMutation, DNAMethylation, CopyNumberVariation

from .proteomics import Protein

from .clinical import ClinicalData

from .multiomics import MultiOmics


def set_backend(new):
    """
    Args:
        new:
    """
    assert new in ["dask", "pandas"]

    if new == "dask":
        this.config["backend"] = dd
    else:
        this.config["backend"] = pd


def set_cache_dir(path, delete_temp=False):
    if not os.path.exists(path):
        raise NotADirectoryError(path)

    this.config["cache_dir"] = path
    astropy.config.set_temp_cache(path=path, delete=delete_temp)
    logging.info("Cache directory is", this.config["cache_dir"])
