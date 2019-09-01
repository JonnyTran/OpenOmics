# -*- coding: utf-8 -*-
"""Tests for `openomics` package."""

import pytest
import os
import random

import openomics as oo
from openomics.multiomics import MultiOmicsData
from openomics.transcriptomics import MessengerRNA, MicroRNA, LncRNA
from openomics.genomics import SomaticMutation

cohort_folder_path = "tests/data/TCGA_LUAD"

@pytest.fixture
def generate_TCGA_LUAD_MessengerRNA():
    return MessengerRNA("LUAD", file_path=os.path.join(cohort_folder_path, "LUAD__geneExp.txt"),
                        index="gene_name", columns="GeneSymbol|TCGA", genes_col_name="GeneSymbol")

@pytest.fixture
def generate_TCGA_LUAD_MicroRNA():
    return MicroRNA("LUAD", file_path=os.path.join(cohort_folder_path, "LUAD__miRNAExp__RPM.txt"),
                    index="gene_name", columns="GeneSymbol|TCGA", genes_col_name="GeneSymbol")

@pytest.fixture
def generate_TCGA_LUAD_LncRNA():
    return LncRNA("LUAD", file_path=os.path.join(cohort_folder_path, "TCGA-rnaexpr.tsv"),
                  index="gene_id", columns="Gene_ID|TCGA", genes_col_name="Gene_ID")

@pytest.fixture
def generate_TCGA_LUAD_SomaticMutation():
    return SomaticMutation("LUAD", file_path=os.path.join(cohort_folder_path, "LUAD__somaticMutation_geneLevel.txt"),
                           index="gene_id", columns="GeneSymbol|TCGA", genes_col_name="GeneSymbol")

def test_import_expression_table_size(generate_TCGA_LUAD_MessengerRNA):
    cohort_name = "LUAD"
    luad_data = MultiOmicsData(cohort_name, import_clinical=True, clinical_file=os.path.join(cohort_folder_path, "nationwidechildrens.org_clinical_patient_luad.txt"))
    luad_data.add_omic(generate_TCGA_LUAD_MessengerRNA)
    luad_data.build_samples()
    assert luad_data.data["MessengerRNA"].shape == generate_TCGA_LUAD_MessengerRNA.expressions.shape

@pytest.fixture
def generate_TCGA_LUAD(generate_TCGA_LUAD_MessengerRNA, generate_TCGA_LUAD_MicroRNA, generate_TCGA_LUAD_LncRNA):
    cohort_name = "LUAD"
    luad_data = MultiOmicsData(cohort_name, import_clinical=True, clinical_file=os.path.join(cohort_folder_path, "nationwidechildrens.org_clinical_patient_luad.txt"))
    luad_data.add_omic(generate_TCGA_LUAD_MessengerRNA)
    luad_data.add_omic(generate_TCGA_LUAD_MicroRNA)
    luad_data.add_omic(generate_TCGA_LUAD_LncRNA)
    return luad_data

def test_TCGA_LUAD_multiomics_transcriptomics(generate_TCGA_LUAD):
    assert all(elem in generate_TCGA_LUAD.get_omics_list()  for elem in [MessengerRNA.name(), MicroRNA.name(), LncRNA.name()])
