import pandas as pd
import os

from TCGAMultiOmics.clinical import ClinicalData
from TCGAMultiOmics.genomic import GeneExpression, SomaticMutation, DNAMethylation, MiRNAExpression, CopyNumberVariation, \
    ProteinExpression, LncRNAExpression


# from TCGAMultiOmics.slideimage import WholeSlideImages


class MultiOmicsData:

    def __init__(self, cancer_type, folder_path, modalities=["WSI", "GE", "SNP", "CNV", "DNA", "MIR", "PRO"]):
        """
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
                    TargetScan/
                        miR_Family_Info.txt
                        Predicted_Targets_Context_Scores.default_predictions.txt
                        Predicted_Targets_Info.default_predictions.txt
                cnv/
                    copyNumber.txt
                protein_rppa/
                    protein_RPPA.txt
                somatic/
                    somaticMutation_geneLevel.txt
                lncrna/
                    TCGA-rnaexpr.tsv
                    HGNC_RNA_long_non-coding.txt

        :param cancer_type: TCGA cancer code name
        :param folder_path: relative directory path to the folder containing multi-omics data downloaded from TCGA-assembler
        """
        self.cancer_type = cancer_type
        self.modalities = modalities

        # LOADING DATA FROM FILES
        self.data = {}
        self.clinical = ClinicalData(cancer_type, os.path.join(folder_path, "clinical/"))
        self.data["PATIENTS"] = self.clinical.patient
        # self.multi_omics_data["BIOSPECIMENS"] = self.clinical.biospecimen
        self.data["DRUGS"] = self.clinical.drugs

        if ("WSI" in modalities):
            pass
            # self.WSI = WholeSlideImages(cancer_type, folder_path)
            # self.data["WSI"] = self.WSI
        if ("GE" in modalities):
            self.GE = GeneExpression(cancer_type, os.path.join(folder_path, "gene_exp/"))
            self.data["GE"] = self.GE.data
        if ("SNP" in modalities):
            self.SNP = SomaticMutation(cancer_type, os.path.join(folder_path, "somatic/"))
            self.data["SNP"] = self.SNP.data
        if ("MIR" in modalities):
            self.MIR = MiRNAExpression(cancer_type, os.path.join(folder_path, "mirna/"))
            self.data["MIR"] = self.MIR.data

            try:
                self.MIR.process_target_scan(mirna_list=self.MIR.get_genes_list(),
                                             gene_symbols=self.GE.get_genes_list(),
                                             targetScan_miR_family_info_path=os.path.join(folder_path, "mirna", "TargetScan", "miR_Family_Info.txt"),
                                             targetScan_predicted_targets_path=os.path.join(folder_path, "mirna", "TargetScan", "Predicted_Targets_Info.default_predictions.txt"),
                                             targetScan_predicted_targets_context_score_path=os.path.join(folder_path, "mirna", "TargetScan", "Predicted_Targets_Context_Scores.default_predictions.txt"))
            except Exception as e:
                print(e)
                print("Could not run MiRNAExpression.process_target_scan() because of missing TargetScan files in directory mirna/TargetScan/", folder_path)

        if ("LNC" in modalities):
            self.LNC = LncRNAExpression(cancer_type, os.path.join(folder_path, "lncrna/"),
                                        HGNC_lncRNA_names_path=os.path.join(folder_path, "lncrna/", "HGNC_RNA_long_non-coding.txt"))
            self.data["LNC"] = self.LNC.data
        if ("DNA" in modalities):
            self.DNA = DNAMethylation(cancer_type, os.path.join(folder_path, "dna/"))
            self.data["DNA"] = self.DNA.data
        if ("CNV" in modalities):
            self.CNV = CopyNumberVariation(cancer_type, os.path.join(folder_path, "cnv/"))
            self.data["CNV"] = self.CNV.data
        if ("PRO" in modalities):
            self.PRO = ProteinExpression(cancer_type, os.path.join(folder_path, "protein_rppa/"))
            self.data["PRO"] = self.PRO.data

        # Build a table for each samples's clinical data
        all_samples = pd.Index([])
        for modality in self.modalities:
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
        Query and fetch the multi-omics dataset based on requested . The data matrices are row index-ed by sample barcode.

        :param modalities: A list of the data modalities to load. Default "all" to select all modalities
        :param target: The clinical data field to include in the
        :param pathologic_stages: List. Only fetch samples having certain stages in their corresponding patient's clinical
        data. For instance, ["Stage I", "Stage II"] will only fetch samples from Stage I and Stage II patients. Default is [] which fetches all pathologic stages.
        :param histological_subtypes: A list specifying the histological subtypes to fetch. Default is [] which fetches all histological sybtypes.
        :param predicted_subtypes: A list specifying the predicted subtypes (if not null) to fetch. Default is [] which fetches all predicted subtypes.
        :param tumor_normal: ["Tumor"] or ["Normal"]. Default is [], which fetches all tumor or normal sample types.
        :param samples_barcode: A list of sample's barcode. If not None, only fetch data with matching bcr_sample_barcodes provided in this list
        :return: X, y
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
            X_multiomics[modality] = self.data[modality].loc[matched_samples, :]

        return X_multiomics, y

    def get_patients_clinical(self, matched_samples):
        """
        Fetch patient's clinical data for each given samples barcodes in the matched_samples
        :param matched_samples: A list of sample barcodes
        """
        return self.data["SAMPLES"].reindex(matched_samples)


    def print_sample_sizes(self):
        for modality in self.data.keys():
            print(modality, self.data[modality].shape if hasattr(self.data[modality],
                                                                             'shape') else "Didn't import data")


    def add_subtypes_to_patients_clinical(self, dictionary):
        """
        This function adds a "predicted_subtype" field to the patients clinical data. For instance, patients were classified
        into subtypes based on their expression profile using k-means, then, to use this function, do:

        add_subtypes_to_patients_clinical(dict(zip(<list of patient barcodes>, <list of corresponding patient's subtypes>)))

        Adding a field to the patients clinical data allows TCGAMultiOmics to query the patients data through the
        .load_data(predicted_subtypes=[]) parameter,

        :param dictionary: A dictionary mapping patient's barcode to a subtype
        """
        self.data["PATIENTS"] = self.data["PATIENTS"].assign(
            predicted_subtype=self.data["PATIENTS"]["bcr_patient_barcode"].map(dictionary))

    def add_association(self, modality_A, modality_B, bi_direction=True):
        pass

