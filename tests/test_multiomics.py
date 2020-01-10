# -*- coding: utf-8 -*-
"""Tests for `openomics` package."""

import os

import pytest

from openomics import MessengerRNA, MicroRNA, LncRNA
from openomics import MultiOmics
from openomics import SomaticMutation

cohort_folder_path = "tests/data/TCGA_LUAD"

@pytest.fixture
def generate_TCGA_LUAD_MessengerRNA():
    return MessengerRNA("LUAD", file_path=os.path.join(cohort_folder_path, "LUAD__geneExp.txt"),
                        columns="GeneSymbol|TCGA", genes_col_name="GeneSymbol", gene_index_by="gene_name",
                        sample_index_by=None)

@pytest.fixture
def generate_TCGA_LUAD_MicroRNA():
    return MicroRNA("LUAD", file_path=os.path.join(cohort_folder_path, "LUAD__miRNAExp__RPM.txt"),
                    columns="GeneSymbol|TCGA", genes_col_name="GeneSymbol", gene_index_by="gene_name",
                    sample_index_by=None)

@pytest.fixture
def generate_TCGA_LUAD_LncRNA():
    return LncRNA("LUAD", file_path=os.path.join(cohort_folder_path, "TCGA-rnaexpr.tsv"), columns="Gene_ID|TCGA",
                  genes_col_name="Gene_ID", gene_index_by="gene_id", sample_index_by=None)

@pytest.fixture
def generate_TCGA_LUAD_SomaticMutation():
    return SomaticMutation("LUAD", file_path=os.path.join(cohort_folder_path, "LUAD__somaticMutation_geneLevel.txt"),
                           columns="GeneSymbol|TCGA", genes_col_name="GeneSymbol", gene_index_by=None,
                           sample_index_by=None)

def test_import_expression_table_size(generate_TCGA_LUAD_MessengerRNA):
    cohort_name = "LUAD"
    luad_data = MultiOmics(cohort_name)
    luad_data.add_clinical_data(
        clinical_data=os.path.join(cohort_folder_path, "nationwidechildrens.org_clinical_patient_luad.txt"))
    luad_data.add_omic(generate_TCGA_LUAD_MessengerRNA)
    luad_data.build_samples()
    print(luad_data.data.keys())
    assert luad_data.data[MessengerRNA.name()].shape == generate_TCGA_LUAD_MessengerRNA.expressions.shape

@pytest.fixture
def generate_TCGA_LUAD(generate_TCGA_LUAD_MessengerRNA, generate_TCGA_LUAD_MicroRNA, generate_TCGA_LUAD_LncRNA):
    cohort_name = "LUAD"
    luad_data = MultiOmics(cohort_name)
    luad_data.add_clinical_data(
        clinical_data=os.path.join(cohort_folder_path, "nationwidechildrens.org_clinical_patient_luad.txt"))
    luad_data.add_omic(generate_TCGA_LUAD_MessengerRNA)
    luad_data.add_omic(generate_TCGA_LUAD_MicroRNA)
    luad_data.add_omic(generate_TCGA_LUAD_LncRNA)
    luad_data.MicroRNA.expressions.index = luad_data.MicroRNA.expressions.index.str.replace("mir", "miR")
    luad_data.MicroRNA.annotations.index = luad_data.MicroRNA.annotations.index.str.replace("mir", "miR")
    return luad_data

def test_TCGA_LUAD_multiomics_transcriptomics(generate_TCGA_LUAD):
    assert all(elem in generate_TCGA_LUAD.get_omics_list()  for elem in [MessengerRNA.name(), MicroRNA.name(), LncRNA.name()])
