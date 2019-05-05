.. openTCGA documentation master file, created by
   sphinx-quickstart on Sat May  4 22:39:31 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to openTCGA's documentation!
====================================

This Python package provide a series of tool to integrate and query the genomics, transcriptomics, proteomics, and clinical TCGA data.
By providing a series of data manipulation tools, openTCGA facilitates the common coding tasks when preparing data for bioinformatics analysis.

Installation via pip (Python >= 3.6.0):

    pip install git+https://github.com/JonnyTran/openTCGA


The TCGA multi-omics data is downloaded from [TCGA-Assembler](http://www.compgenome.org/TCGA-Assembler/).
Load all multi-omics data files according the following folder structure and naming convention:

    tcga_data_path/
        clinical/
            genome.wustl.edu_biospecimen_sample.txt (optional)
            nationwidechildrens.org_clinical_drug.txt
            nationwidechildrens.org_clinical_patient.txt
        gene_exp/
            geneExp.txt
        mirna/
            miRNAExp__RPM.txt
        lncrna/
            TCGA-rnaexpr.tsv
        cnv/
            copyNumber.txt
        protein_rppa/
            protein_RPPA.txt
        somatic/
            somaticMutation_geneLevel.txt

The microRNA and lncRNA data requires additional external databases, e.g. TargetScan, microRNA family, HGNC long non-coding RNA names, etc.

    external_data_path/
        TargetScan/
            Gene_info.txt
            miR_Family_Info.txt
            Predicted_Targets_Context_Scores.default_predictions.txt
            Predicted_Targets_Info.default_predictions.txt

        HUGO_Gene_names/
            gene_with_protein_product.txt
            RNA_long_non-coding.txt
            RNA_micro.txt


# How to use openTCGA:


## Importing the openTCGA library


```python
from openTCGA.multiomics import MultiOmicsData
```

## Import TCGA LUAD data downloaded from TCGA-Assembler


```python
folder_path ="./data/tcga-assembler/LUAD/"
```



```python
# Load all modalities: Gene Expression, MicroRNA expression lncRNA expression, Copy Number Variation, Somatic Mutation, DNA Methylation, and Protein Expression data
luad_data = MultiOmicsData(cancer_type="LUAD", folder_path=folder_path,
                           modalities=["GE", "MIR", "LNC", "CNV", "SNP", "PRO"])

```

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
