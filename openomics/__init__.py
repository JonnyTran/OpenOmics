# -*- coding: utf-8 -*-

"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = 'nhat.tran@mavs.uta.edu'
__version__ = '0.7.3'

from __future__ import print_function, division, absolute_import

try:
    from .import (
        clinical,
        genomics,
        multiomics,
        multicohorts,
        transcriptomics,
        proteomics,
        database,
        utils
    )

    from .database import (
        annotation,
        disease,
        interaction,
        ontology
    )

    from .utils import (
        io
    )
except ImportError as e:
    msg = (
        "OpenOmics requirements are not installed.\n\n"
        "Please pip install as follows:\n\n"
        "  pip install openomics --upgrade  # or pip install"
    )
    raise ImportError(str(e) + "\n\n" + msg)
