from __future__ import print_function, division, absolute_import

import json
import logging
import os
import sys

import astropy
import dask.dataframe as dd
import pandas as pd

"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = 'nhat.tran@mavs.uta.edu'
__version__ = '0.8.5'

config_file_path = "~/.openomics/conf.json"

# Initialize configurations
this = sys.modules[__name__]
this.config = {}

# Set pandas backend
this.config["backend"] = pd

# Set cache download directory
this.config["cache_dir"] = astropy.config.get_cache_dir(this.__name__)
logging.info("Cache directory is", this.config["cache_dir"])

# Initialize user configuration file at ~/.openomics/conf.json
if not os.path.exists(config_file_path):
    if not os.path.exists("~/.openomics/"):
        os.mkdir("~/.openomics/")

    if not os.path.exists(config_file_path):
        base_config = {}
        base_config['database'].append({
            'cache_dir': astropy.config.get_cache_dir(this.__name__)
        })

        with open(config_file_path, 'w') as config_file:
            json.dump(base_config, config_file)

# Read configuration from ~/.openomics/conf.json
if os.path.exists(config_file_path):
    user_config = json.load(config_file_path)
    if user_config:
        this.config.update(user_config)

from . import database, utils

from .transcriptomics import (
    Expression, MessengerRNA, MicroRNA, LncRNA,
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
