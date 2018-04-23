import pandas as pd

import os

from definitions import ROOT_DIR

from clinicaldata import ClinicalData
from genomicdata import GeneExpression, SNP, DNAMethylation, miRNAExpression, CopyNumberVariation, \
    ProteinExpression, LncRNAExpression


# from src.features.slide_image import WholeSlideImages


class MultiOmicsData:

    def __init__(self, cancer_type, folder_path, modalities=["WSI", "GE", "SNP", "CNV", "DNA", "MIR", "PRO"]):
        """
        Load all multi-omics TCGA data from a given folder_path with the following folder structure:

            folder_path/
                clinical/
                    genome.wustl.edu_biospecimen_sample.txt
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

        :param cancer_type: TCGA cancer code name
        :param folder_path: relative directory path to the folder containing multi-omics data downloaded from TCGA-assembler
        """
        self.cancer_type = cancer_type
        self.modalities = modalities

        # LOADING DATA FROM FILES
        self.data = {}
        self.clinical = ClinicalData(cancer_type, folder_path + "clinical/")
        self.data["PATIENTS"] = self.clinical.patient
        # self.multi_omics_data["BIOSPECIMENS"] = self.clinical.biospecimen
        self.data["DRUGS"] = self.clinical.drugs

        if ("WSI" in modalities):
            self.WSI = WholeSlideImages(cancer_type, folder_path)
            self.data["WSI"] = self.WSI
        if ("GE" in modalities):
            self.GE = GeneExpression(cancer_type, folder_path + "gene_exp/", )
            self.data["GE"] = self.GE.data
        if ("SNP" in modalities):
            self.SNP = SNP(cancer_type, folder_path + "somatic/")
            self.data["SNP"] = self.SNP.data
        if ("MIR" in modalities):
            self.MIR = miRNAExpression(cancer_type, folder_path + "mirna/")
            self.data["MIR"] = self.MIR.data

            self.MIR.process_target_scan(mirna_list=self.MIR.get_genes_list(),
                                         gene_symbols=self.GE.get_genes_list())
        if ("LNC" in modalities):
            self.LNC = LncRNAExpression(cancer_type, folder_path + "lncrna/")
            self.data["LNC"] = self.LNC.data
        if ("DNA" in modalities):
            self.DNA = DNAMethylation(cancer_type, folder_path + "dna/")
            self.data["DNA"] = self.DNA.data
        if ("CNV" in modalities):
            self.CNV = CopyNumberVariation(cancer_type, folder_path + "cnv/")
            self.data["CNV"] = self.CNV.data
        if ("PRO" in modalities):
            self.PRO = ProteinExpression(cancer_type, folder_path + "protein_rppa/")
            self.data["PRO"] = self.PRO.data

        # Build a table for each samples's clinical data
        all_samples = pd.Index([])
        for modality in modalities:
            all_samples = all_samples.union(self.data[modality].index)
        self.clinical.build_clinical_samples(all_samples)
        self.data["SAMPLES"] = self.clinical.samples

        self.print_sample_sizes()

    def match_samples(self, modalities):
        """
        Return the index of bcr_sample_barcodes of the intersection of samples from all modalities

        :param modalities: An array of modalities
        :return: An pandas Index list
        """
        # TODO check that for single modalities, this fetch all patients
        matched_samples = self.data[modalities[0]].index.copy()

        for modality in modalities:
            matched_samples = matched_samples.join(self.data[modality].index, how="inner")

        return matched_samples

    def load_data(self, modalities, target=['pathologic_stage'],
                  pathologic_stages=[], histological_subtypes=[], predicted_subtypes=[], tumor_normal=[],
                  samples_barcode=None):
        """
        Load and return the multi-omics dataset (classification)
        :param modalities: A list of the data modalities to load. Default "all" to select all modalities
        """
        if modalities == 'all' or modalities == None:
            modalities = self.modalities
        elif modalities:
            modalities = modalities
        else:
            raise Exception("Need to specify which multi-omics to retrieve")

        matched_samples = self.match_samples(modalities)
        if not (samples_barcode is None):
            matched_samples = samples_barcode

        # Build targets clinical data
        y = self.get_patients_clinical(matched_samples)

        # Select only samples with certain cancer stage or subtype
        if pathologic_stages:
            y = y[y['pathologic_stage'].isin(pathologic_stages)]
        if histological_subtypes:
            y = y[y['histologic_subtype'].isin(histological_subtypes)]
        if predicted_subtypes:
            y = y[y['predicted_subtype'].isin(predicted_subtypes)]
        if tumor_normal:
            y = y[y['tumor_normal'].isin(tumor_normal)]
        # TODO if normal_matched:
        #     target =

        # Filter y target column labels
        y = y.filter(target)
        y.dropna(axis=0, inplace=True)

        matched_samples = y.index

        # Build data matrix for each modality, indexed by matched_samples
        X_multiomics = {}
        for modality in modalities:
            X_multiomics[modality] = self.data[modality].loc[matched_samples]

        return X_multiomics, y

    def get_patients_clinical(self, matched_samples):
        """
        Fetch patient's clinical data for each given samples barcodes given in matched_samples
        """
        return self.data["SAMPLES"].loc[matched_samples]


    def print_sample_sizes(self):
        for modality in self.data.keys():
            print(modality, self.data[modality].shape if hasattr(self.data[modality],
                                                                             'shape') else "Didn't import data")


    def add_subtypes_to_patients_clinical(self, dictionary):
        self.data["PATIENTS"] = self.data["PATIENTS"].assign(
            predicted_subtype=self.data["PATIENTS"]["bcr_patient_barcode"].map(dictionary))

    def add_association(self, modality_A, modality_B, bi_direction=True):
        pass


if __name__ == '__main__':
    folder_path = "/data/tcga-assembler/LUSC/"
    luad_data = MultiOmicsData(cancer_type="LUSC", folder_path=ROOT_DIR + folder_path,
                               modalities=["GE", "MIR"])

    matched_samples = luad_data.match_samples(modalities=["GE"])
    # print("matched samples", matched_samples.shape, matched_samples)


    patients_clinical = luad_data.get_patients_clinical(matched_samples)
    X, y = luad_data.load_data(modalities="all", target=['pathologic_stage'])
    # print(patients_clinical)
