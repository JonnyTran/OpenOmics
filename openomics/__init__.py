from __future__ import print_function, division, absolute_import
# -*- coding: utf-8 -*-

"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = 'nhat.tran@mavs.uta.edu'
__version__ = '0.7.7'

try:
    from .clinical import ClinicalData

    from .genomics import (
        SomaticMutation, DNAMethylation, CopyNumberVariation
    )

    from .transcriptomics import (
        MessengerRNA, MicroRNA, LncRNA
    )

    from .proteomics import (
        Protein
    )

    from .multiomics import (
        MultiOmics
    )

    # from .database import (
    #     annotation,
    #     disease,
    #     interaction,
    #     ontology
    # )


except ImportError as e:
    msg = (
        "OpenOmics requirements are not installed.\n\n"
        "Please pip install as follows:\n\n"
        "  pip install openomics --upgrade  # or pip install"
    )
    raise ImportError(str(e) + "\n\n" + msg)
