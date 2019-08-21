import os
import pandas as pd

from openTCGA.clinical import ClinicalData, HISTOLOGIC_SUBTYPE, PATHOLOGIC_STAGE, BCR_PATIENT_BARCODE, \
    TUMOR_NORMAL, PREDICTED_SUBTYPE
from openTCGA.expression import GeneExpression, MiRNAExpression, \
    ProteinExpression, LncRNAExpression
from openTCGA.genomic import SomaticMutation, DNAMethylation, CopyNumberVariation
from openTCGA.slideimage import WholeSlideImages

class MultiOmicsData:
    def __init__(self, cohort_name:str, cohort_folder_path:str, external_data_path:str, modalities:list, import_sequences="longest",
                 replace_U2T=True, remove_duplicate_genes=True, auto_import_clinical=True, process_genes_info=True):
        """
        .. class:: MultiOmicsData
        Load all multi-omics TCGA data from a given tcga_data_path with the following folder structure:
            cohort_folder/
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
                wsi/
                    ...

        Load the external data downloaded from various databases. These data will be imported as attribute information to
        the genes, or interactions between the genes.

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

        :param cohort_name: TCGA cancer cohort name
        :param cohort_folder_path: directory path to the folder containing clinical and multi-omics data downloaded from TCGA-assembler
        :param external_data_path: directory path to the folder containing external databases
        :param modalities: A list of multi-omics data to import. All available data includes ["CLI", "WSI", "GE", "SNP", "CNV", "DNA", "MIR", "LNC", "PRO"]. Clinical data is always automatically imported.
        """
        self.cancer_type = cohort_name
        self.modalities = modalities
        self.external_data_path = external_data_path

        # LOADING DATA FROM FILES
        self.data = {}
        if auto_import_clinical or ("CLI" in modalities):
            self.clinical = ClinicalData(cohort_name, os.path.join(cohort_folder_path, "clinical/"))
            self.data["PATIENTS"] = self.clinical.patient
            if hasattr(self.clinical, "biospecimen"):
                self.data["BIOSPECIMENS"] = self.clinical.biospecimen
            if hasattr(self.clinical, "drugs"):
                self.data["DRUGS"] = self.clinical.drugs

        if "WSI" in modalities:
            self.WSI = WholeSlideImages(cohort_name, os.path.join(cohort_folder_path, "wsi/"))
            self.data["WSI"] = self.WSI

        if "GE" in modalities:
            GE_file_path = os.path.join(cohort_folder_path, "gene_exp", "geneExp.txt")
            self.GE = GeneExpression(cohort_name, GE_file_path,
                                     import_sequences=import_sequences, replace_U2T=replace_U2T)
            self.data["GE"] = self.GE.expression

            try:
                self.GE.process_targetScan_gene_info(
                    targetScan_gene_info_path=os.path.join(external_data_path, "TargetScan", "Gene_info.txt"))

                self.GE.process_HUGO_protein_coding_genes_info(
                    os.path.join(external_data_path, "HUGO_Gene_names", "gene_with_protein_product.txt"))

                self.GE.process_GO_genes_info(os.path.join(external_data_path, "GeneOntology"))

                self.GE.process_genemania_interactions(os.path.join(external_data_path, "GeneMania"))

                self.GE.process_GENCODE_transcript_data(gencode_folder_path=os.path.join(external_data_path, "GENCODE"))

                self.GE.process_biogrid_GRN_edgelist(biogrid_folder_path=os.path.join(external_data_path, "BioGRID"))

                self.GE.process_RegNet_gene_regulatory_network(
                    grn_file_path=os.path.join(external_data_path, "RegNetwork", "human.source"))

                self.GE.process_DisGeNET_gene_disease_associations(
                    disgenet_folder_path=os.path.join(external_data_path, "DisGeNET"))

                self.GE.process_starBase_RNA_RNA_interactions(os.path.join(external_data_path, "StarBase v2.0"))
            except FileNotFoundError as e:
                print(e)

        if "SNP" in modalities:
            file_path_SNP = os.path.join(cohort_folder_path, "somatic/", "somaticMutation_geneLevel.txt")
            self.SNP = SomaticMutation(cohort_name, file_path_SNP)
            self.data["SNP"] = self.SNP.expression

        if "MIR" in modalities:
            file_path_MIR = os.path.join(cohort_folder_path, "mirna/", "miRNAExp__RPM.txt")
            self.MIR = MiRNAExpression(cohort_name, file_path_MIR,
                                       import_sequences=import_sequences, replace_U2T=replace_U2T)
            self.data["MIR"] = self.MIR.expression

            try:
                self.MIR.process_mirbase_data(mirbase_folder_path=os.path.join(external_data_path, "mirbase"))
                self.MIR.process_target_scan(targetScan_folder_path=os.path.join(external_data_path, "TargetScan"))

                self.MIR.process_miRTarBase_miRNA_target_interactions(
                    miRTarBase_path=os.path.join(external_data_path, "miRTarBase"))

                self.MIR.process_mirnadisease_associations(
                    HMDD_miRNAdisease_path=os.path.join(external_data_path, "HMDD_miRNAdisease"))

                self.MIR.process_HUGO_miRNA_gene_info(
                    HUGO_folder_path=os.path.join(external_data_path, "HUGO_Gene_names"))
                self.MIR.process_RNAcentral_annotation_info(
                    RNAcentral_folder_path=os.path.join(external_data_path, "RNAcentral"))

            except FileNotFoundError as e:
                print(e)
                print(
                    "Could not run MiRNAExpression.process_target_scan() because of missing TargetScan data folder in the directory",
                    external_data_path)

        if "LNC" in modalities:
            file_path_LNC = os.path.join(cohort_folder_path, "lncrna/", "TCGA-rnaexpr.tsv")
            self.LNC = LncRNAExpression(cohort_name, file_path_LNC,
                                        HGNC_lncRNA_names_file_path=os.path.join(external_data_path, "HUGO_Gene_names",
                                                                                 "RNA_long_non-coding.txt"),
                                        GENCODE_folder_path=os.path.join(external_data_path, "GENCODE"),
                                        external_data_path=external_data_path, import_sequences=import_sequences,
                                        replace_U2T=replace_U2T)
            self.data["LNC"] = self.LNC.expression

            try:
                self.LNC.process_lncRNome_miRNA_binding_sites(os.path.join(external_data_path, "lncRNome"))
                self.LNC.process_lncRNome_gene_info(os.path.join(external_data_path, "lncRNome"))
                self.LNC.process_lncBase_miRNA_lncRNA_interactions(
                    lncBase_folder_path=os.path.join(external_data_path, "lncBase"))
                self.LNC.process_starBase_miRNA_lncRNA_interactions(os.path.join(external_data_path, "StarBase v2.0"))
                self.LNC.process_starBase_lncRNA_RNA_interactions(os.path.join(external_data_path, "StarBase v2.0"))
                self.LNC.process_LncReg_lncRNA_RNA_regulatory_interactions(
                    LncReg_folder_path=os.path.join(external_data_path, "LncReg"))
                self.LNC.process_lncrna2target_interactions(os.path.join(external_data_path, "lncrna2target"))
                self.LNC.process_lncRInter_interactions(os.path.join(external_data_path, "lncRInter"))
                self.LNC.process_NPInter_ncRNA_RNA_regulatory_interactions(
                    NPInter_folder_path=os.path.join(external_data_path, "NPInter"))
                self.LNC.process_NONCODE_func_annotation(os.path.join(external_data_path, "NONCODE"))
                self.LNC.process_lncrnadisease_associations(
                    lncrnadisease_folder_path=os.path.join(external_data_path, "lncrnadisease"))
                self.LNC.process_RNAcentral_annotation_info(
                    RNAcentral_folder_path=os.path.join(external_data_path, "RNAcentral"))
            except FileNotFoundError as e:
                print(e)

        if "DNA" in modalities:
            file_path_DNA = os.path.join(cohort_folder_path, "dna/", "methylation_450.txt")
            self.DNA = DNAMethylation(cohort_name, file_path_DNA)
            self.data["DNA"] = self.DNA.expression

        if "CNV" in modalities:
            file_path_CNV = os.path.join(cohort_folder_path, "cnv/", "copyNumber.txt")
            self.CNV = CopyNumberVariation(cohort_name, file_path_CNV)
            self.data["CNV"] = self.CNV.expression

        if "PRO" in modalities:
            file_path_PRO = os.path.join(cohort_folder_path, "protein_rppa/", "protein_RPPA.txt")
            self.PRO = ProteinExpression(cohort_name, file_path_PRO)
            self.data["PRO"] = self.PRO.expression
            self.PRO.process_HPRD_PPI_network(
                ppi_data_file_path=os.path.join(external_data_path, "HPRD_PPI",
                                                "BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt"))

        # Build a table for each samples's clinical data
        if len(modalities) > 1:  # TODO Has to make sure at least one GenomicData present
            all_samples = pd.Index([])
            for modality in self.modalities:
                all_samples = all_samples.union(self.data[modality].index)
            self.clinical.build_clinical_samples(all_samples)
            self.data["SAMPLES"] = self.clinical.samples

        # Remove duplicate genes between different multi-omics (e.g. between gene expression and lncRNA expressions
        if remove_duplicate_genes:
            if "GE" in modalities and "LNC" in modalities:
                self.GE.drop_genes(set(self.GE.get_genes_list()) & set(self.LNC.get_genes_list()))

        self.print_sample_sizes()

        if process_genes_info:
            for modality in modalities:
                if hasattr(self[modality], "process_genes_info"):
                    self[modality].process_genes_info()
                    print("Processed genes info for ", modality)

    def __getitem__(self, item):
        """
        This function allows the MultiOmicData class objects to access individual omics by a dictionary lookup
        """
        if item == "GE":
            return self.GE
        elif item == "MIR":
            return self.MIR
        elif item == "LNC":
            return self.LNC
        elif item == "WSI":
            return self.WSI
        elif item == "SNP":
            return self.SNP
        elif item == "CNV":
            return self.CNV
        elif item == "DNA":
            return self.DNA
        elif item == "PRO":
            return self.PRO
        elif item == "CLI":
            return self.clinical.patient
        elif item == "DRU":
            return self.clinical.drugs

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

    # noinspection PyPep8Naming
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
        if modalities == 'all' or modalities is None:
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
            y = y[y[PATHOLOGIC_STAGE].isin(pathologic_stages)]
        if histological_subtypes:
            y = y[y[HISTOLOGIC_SUBTYPE].isin(histological_subtypes)]
        if predicted_subtypes:
            y = y[y[PREDICTED_SUBTYPE].isin(predicted_subtypes)]
        if tumor_normal:
            y = y[y[TUMOR_NORMAL].isin(tumor_normal)]
        # TODO if normal_matched:
        #     target =

        # Filter y target column labels
        y = y.filter(target)
        y.dropna(axis=0, inplace=True)

        matched_samples = y.index

        # Build data matrix for each modality, indexed by matched_samples
        X_multiomics = {}
        for modality in modalities:
            X_multiomics[modality] = self.data[modality].loc[matched_samples, self[modality].get_genes_list()]

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

        Adding a field to the patients clinical data allows openTCGA to query the patients data through the
        .load_data(predicted_subtypes=[]) parameter,

        :param dictionary: A dictionary mapping patient's barcode to a subtype
        """
        self.data["PATIENTS"] = self.data["PATIENTS"].assign(
            predicted_subtype=self.data["PATIENTS"][BCR_PATIENT_BARCODE].map(dictionary))
