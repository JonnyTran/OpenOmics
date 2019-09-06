from typing import List, Dict, Union

import pandas as pd

from openomics.clinical import ClinicalData, HISTOLOGIC_SUBTYPE, PATHOLOGIC_STAGE, TUMOR_NORMAL, PREDICTED_SUBTYPE
from openomics.genomics import SomaticMutation, CopyNumberVariation, DNAMethylation
from openomics.image import WholeSlideImage
from openomics.proteomics import Protein
from openomics.transcriptomics import MessengerRNA, MicroRNA, LncRNA, ExpressionData


class MultiOmics:
    def __init__(self, cohort_name, import_clinical=True, clinical_file=None):
        # type: (str, List[str], bool, str) -> None
        """
        Load all multi-omics data from a given cohort_folder path.


        Args:
            cohort_name (str): the clinical cohort name
            import_clinical (bool, ClinicalData):
        """
        self.cancer_type = cohort_name
        self.omics_list = []

        # This is a data dictionary accessor to retrieve DataFrame's
        self.data = {}

        if import_clinical and type(import_clinical) == ClinicalData:
            self.clinical = import_clinical
        elif import_clinical:
            self.clinical = ClinicalData(cohort_name, clinical_file)

        if import_clinical:
            self.data["PATIENTS"] = self.clinical.patient
            if hasattr(self.clinical, "biospecimen"):
                self.data["BIOSPECIMENS"] = self.clinical.biospecimen
            if hasattr(self.clinical, "drugs"):
                self.data["DRUGS"] = self.clinical.drugs

        # if "GE" in omics:
        #     table_path_GE = os.path.join(cohort_folder, "gene_exp", "geneExp.txt")
        #     self.GE = MessengerRNA(cohort_name, table_path_GE, columns="GeneSymbol|TCGA", index="GeneSymbol")
        #     self.data["GE"] = self.GE.expressions

            # try:
            #
            #     self.GE.process_HUGO_protein_coding_genes_info(
            #         os.path.join(external_data_path, "HUGO_Gene_names", "gene_with_protein_product.txt"))
            #
            #     self.GE.process_GO_genes_info(os.path.join(external_data_path, "GeneOntology"))
            #
            #     self.GE.process_RegNet_gene_regulatory_network(
            #         grn_file_path=os.path.join(external_data_path, "RegNetwork", "human.source"))
            #
            #     self.GE.process_DisGeNET_gene_disease_associations(
            #         disgenet_folder_path=os.path.join(external_data_path, "DisGeNET"))
            #
            #     self.GE.process_starBase_RNA_RNA_interactions(os.path.join(external_data_path, "StarBase v2.0"))
            # except FileNotFoundError as e:
            #     print(e)

        # if "SNP" in omics:
        #     file_path_SNP = os.path.join(cohort_folder, "somatic/", "somaticMutation_geneLevel.txt")
        #     self.SNP = SomaticMutation(cohort_name, file_path_SNP)
        #     self.data["SNP"] = self.SNP.expressions
        #
        # if "MIR" in omics:
        #     file_path_MIR = os.path.join(cohort_folder, "mirna/", "miRNAExp__RPM.txt")
        #     self.MIR = MicroRNA(cohort_name, file_path_MIR)
        #     self.data["MIR"] = self.MIR.expressions

            # try:
            #     self.MIR.process_mirnadisease_associations(
            #         HMDD_miRNAdisease_path=os.path.join(external_data_path, "HMDD_miRNAdisease"))
            #
            #     # self.MIR.process_HUGO_miRNA_gene_info(
            #     #     HUGO_folder_path=os.path.join(external_data_path, "HUGO_Gene_names"))
            #     # self.MIR.process_RNAcentral_annotation_info(
            #     #     RNAcentral_folder_path=os.path.join(external_data_path, "RNAcentral"))
            #

        # if "LNC" in omics:
        #     file_path_LNC = os.path.join(cohort_folder, "lncrna", "TCGA-rnaexpr.tsv")
        #     self.LNC = LncRNA(cohort_name, file_path_LNC, columns="Gene_ID|TCGA", index="Gene_ID")
        #     self.data["LNC"] = self.LNC.expressions

            # try:
            #     self.LNC.process_lncRNome_miRNA_binding_sites(os.path.join(external_data_path, "lncRNome"))
            #     self.LNC.process_lncRNome_gene_info(os.path.join(external_data_path, "lncRNome"))
            #     self.LNC.process_starBase_miRNA_lncRNA_interactions(os.path.join(external_data_path, "StarBase v2.0"))
            #     self.LNC.process_starBase_lncRNA_RNA_interactions(os.path.join(external_data_path, "StarBase v2.0"))
            #     self.LNC.process_LncReg_lncRNA_RNA_regulatory_interactions(
            #         LncReg_folder_path=os.path.join(external_data_path, "LncReg"))
            #     self.LNC.process_NPInter_ncRNA_RNA_regulatory_interactions(
            #         NPInter_folder_path=os.path.join(external_data_path, "NPInter"))
            #     self.LNC.process_lncrnadisease_associations(
            #         lncrnadisease_folder_path=os.path.join(external_data_path, "lncrnadisease"))

        # if "PRO" in omics:
            # self.PRO.process_HPRD_PPI_network(
            #     ppi_data_file_path=os.path.join(external_data_path, "HPRD_PPI",
            #                                     "BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt"))

        # Build a table for each samples's clinical data
        if import_clinical:
            self.build_samples()

        # Remove duplicate genes between different multi-omics (e.g. between gene expression and lncRNA expressions
        self.print_sample_sizes()

    def remote_duplate_genes(self):
        """
        Removes duplicate genes between any omics such that the index across all omics has no duplicates.
        """
        for omic_A in self.omics_list:
            for omic_B in self.omics_list:
                if omic_A != omic_B:
                    self.__getattribute__(omic_A).drop_genes(set(self.__getattribute__(omic_A).get_genes_list()) & set(
                        self.__getattribute__(omic_B).get_genes_list()))

    def add_omic(self, omic_data, initialize_annotations=True):
        # type: (ExpressionData, bool) -> None
        """
        Adds an omic object to the Multiomics such that the samples in omic matches the samples existing in the other omics.

        Args:
            omic_data (openomics.transcriptomics.ExpressionData):
                The omic to add, e.g., MessengerRNA, MicroRNA, LncRNA, etc.
            initialize_annotations (bool): default True.
                If true, initializes the annotations dataframe in the omic object
        """
        self.__dict__[omic_data.name()] = omic_data

        if omic_data.name not in self.omics_list:
            self.omics_list.append(omic_data.name())

        # dictionary as data accessor to the expression data
        self.data[omic_data.name()] = omic_data.expressions

        # Initialize annotations
        if initialize_annotations:
            omic_data.initialize_annotations(None, omic_data.gene_index)

        print(omic_data.name(), self.data[omic_data.name()].shape if hasattr(self.data[omic_data.name()], 'shape') else ": None")

    def get_omics_list(self):
        return self.omics_list

    def build_samples(self, agg_by="union"):
        """

        Args:
            agg_by (str): ["union", "intersection"]
        """
        if len(self.omics_list) < 1:  # make sure at least one ExpressionData present
            raise Exception("Must add at least one omic to this MultiOmics object.")

        all_samples = pd.Index([])
        for omic in self.omics_list:
            if agg_by == "union":
                all_samples = all_samples.union(self.data[omic].index)
            elif agg_by == "intersection":
                all_samples = all_samples.intersection(self.data[omic].index)

        if hasattr(self, "clinical"):
            self.clinical.build_clinical_samples(all_samples)
            self.data["SAMPLES"] = self.clinical.samples.index
        else:
            self.data["SAMPLES"] = all_samples


    def __getitem__(self, item):
        # type: (str) -> object
        """
        This function allows the MultiOmicData class objects to access individual omics by a dictionary lookup, e.g.

        Args:
            item (str): a string of the class name
        """
        if item.lower() == MessengerRNA.name().lower():
            return self.__getattribute__(MessengerRNA.name())

        elif item.lower() == MicroRNA.name().lower():
            return self.__getattribute__(MicroRNA.name())

        elif item.lower() == LncRNA.name().lower():
            return self.__getattribute__(LncRNA.name())

        elif item.lower() == WholeSlideImage.name().lower():
            return self.__getattribute__(WholeSlideImage.name())

        elif item.lower() == SomaticMutation.name().lower():
            return self.__getattribute__(SomaticMutation.name())

        elif item.lower() == CopyNumberVariation.name().lower():
            return self.__getattribute__(CopyNumberVariation.name())

        elif item.lower() == DNAMethylation.name().lower():
            return self.__getattribute__(DNAMethylation.name())

        elif item.lower() == Protein.name().lower():
            return self.__getattribute__(Protein.name())

        elif item.lower() == ClinicalData.name().lower():
            return self.clinical.patient
        elif item.lower() == "DRU":
            return self.clinical.drugs
        else:
            raise Exception('String accessor must be one of {"MessengerRNA", "MicroRNA", "LncRNA", "Protein", etc.}')

    def match_samples(self, omics):
        """
        Return the index of bcr_sample_barcodes of the intersection of samples from all modalities

        Args:
            omics: An array of modalities
        Returns:
            matched_sapmles: An pandas Index list
        """
        # TODO check that for single modalities, this fetch all patients
        matched_samples = self.data[omics[0]].index.copy()

        for omic in omics:
            matched_samples = matched_samples.join(self.data[omic].index, how="inner")

        return matched_samples

    def load_data(self, omics, target=['pathologic_stage'],
                  pathologic_stages=None, histological_subtypes=None, predicted_subtypes=None, tumor_normal=None,
                  samples_barcode=None):
        # type: (Union[List[str], str], List[str], List[str], List[str], List[str], List[str], List[str]) -> (Dict[str, pd.DataFrame], pd.DataFrame)
        """

        Args:
            omics (list):
                A list of the data modalities to load. Default "all" to select all modalities
            target (list):
                The clinical data fields to include in the
            pathologic_stages (list):
                Only fetch samples having certain stages in their corresponding patient's clinical data. For instance, ["Stage I", "Stage II"] will only fetch samples from Stage I and Stage II patients. Default is [] which fetches all pathologic stages.
            histological_subtypes:
                A list specifying the histological subtypes to fetch. Default is [] which fetches all histological sybtypes.
            predicted_subtypes:
                A list specifying the histological subtypes to fetch. Default is [] which fetches all histological sybtypes.
            tumor_normal: ["Tumor", "Normal"].
                Default is [], which fetches all tumor or normal sample types.
            samples_barcode:
                A list of sample's barcode. If not None, only fetch data with matching samples provided in this list.

        Returns:
            (X, y): Returns X, a dictionary containing the multiomics data that have data

        """

        if omics == 'all' or omics is None:
            omics = self.omics_list
        elif omics:
            omics = omics
        else:
            raise Exception("Need to specify which multi-omics to retrieve")

        matched_samples = self.match_samples(omics)
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

        # Build data matrix for each omic, indexed by matched_samples
        X_multiomics = {}
        for omic in omics:
            X_multiomics[omic] = self.data[omic].loc[matched_samples, self[omic].get_genes_list()]

        return X_multiomics, y

    def get_patients_clinical(self, matched_samples):
        """
        Fetch patient's clinical data for each given samples barcodes in the matched_samples

        Args:
            matched_samples: A list of sample barcodes
        Returns
            samples_index: Index of samples
        """
        return self.data["SAMPLES"].reindex(matched_samples)

    def print_sample_sizes(self):
        for omic in self.data.keys():
            print(omic, self.data[omic].shape if hasattr(self.data[omic],
                                                                 'shape') else "Didn't import data")

    def add_subtypes_to_patients_clinical(self, dictionary):
        """
        This function adds a "predicted_subtype" field to the patients clinical data. For instance, patients were classified
        into subtypes based on their expression profile using k-means, then, to use this function, do:

        add_subtypes_to_patients_clinical(dict(zip(<list of patient barcodes>, <list of corresponding patient's subtypes>)))

        Adding a field to the patients clinical data allows openomics to query the patients data through the
        .load_data(predicted_subtypes=[]) parameter,

        Args:
            dictionary: A dictionary mapping patient's barcode to a subtype
        """
        self.data["PATIENTS"] = self.data["PATIENTS"].assign(
            predicted_subtype=self.data["PATIENTS"][self.clinical.patient_column].map(dictionary))
