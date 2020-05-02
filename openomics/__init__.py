from __future__ import print_function, division, absolute_import

# -*- coding: utf-8 -*-

"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = 'nhat.tran@mavs.uta.edu'
__version__ = '0.7.9'

# from . import database
from .clinical import ClinicalData

from .multiomics import (
    MultiOmics
)

from .visualization import (
    umap
)

from .transcriptomics import (
    ExpressionData, MessengerRNA, MicroRNA, LncRNA,
)

from .genomics import (
    SomaticMutation, DNAMethylation, CopyNumberVariation
)

from .proteomics import (
    Protein
)
