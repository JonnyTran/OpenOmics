# TCGAMultiOmics

This Python package provide a series of tool to download, integrate, and query multi-omics TCGA data. The TCGA data is downloaded from TCGA-Assembler.

Load all multi-omics TCGA data from a given folder_path with the following folder structure:

    folder_path/
        clinical/
            genome.wustl.edu_biospecimen_sample.txt (optional)
            nationwidechildrens.org_clinical_drug.txt
            nationwidechildrens.org_clinical_patient.txt
        gene_exp/
            geneExp.txt
        mirna/
            miRNAExp__RPM.txt
        cnv/
            copyNumber.txt
        protein_rppa/
            protein_RPPA.txt
        somatic/
            somaticMutation_geneLevel.txt
        lncrna/
            TCGA-rnaexpr.tsv