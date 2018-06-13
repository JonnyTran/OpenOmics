import os

import numpy as np
import pandas as pd
import networkx as nx
from TCGAMultiOmics.utils import GTF

class GenomicData:
    def __init__(self, cancer_type, file_path, columns="GeneSymbol|TCGA",
                 import_from_TCGA_Assembler=True, log2_transform=False):
        """

        :param cancer_type: TCGA cancer cohort code name
        :param file_path: Path of the table file to import
        :param columns: column names to import from the table. Columns names imported are string match, separated by "|"
        :param import_from_TCGA_Assembler: If True, perform preprocessing steps for the table data obtained from TCGA-Assembler tool. If False, import a pandas table as-is with bcr_sample_barcode for row index, and gene names as columns
        :param log2_transform: Whether to log2 transform the expression values
        """
        self.cancer_type = cancer_type

        if import_from_TCGA_Assembler:
            self.data = self.preprocess_expression_table(pd.read_table(file_path), columns)

        if log2_transform:
            self.data = self.data.applymap(self.log2_transform)

        # Save samples and features for this omics data
        self.samples = self.data.index
        self.features = self.data.columns.tolist()
        # self.features.remove("bcr_sample_barcode")


    def preprocess_expression_table(self, df, columns):
        """
        This function preprocesses the table files obtained from TCGA-Assembler
        :param df:
        :param columns:
        :return:
        """
        table = df

        # Filter columns
        table = table.filter(regex=columns)

        # Cut TCGA column names to sample barcode
        table.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x, inplace=True)

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        # Drop NA GeneSymbol rows
        table.dropna(axis=0, inplace=True)

        # Remove entries with unknown Gene Symbol
        table = table[table.GeneSymbol != '?']

        # Transpose dataframe to patient rows and GeneSymbol columns
        table.index = table.GeneSymbol
        table.drop(['GeneSymbol'], axis=1, inplace=True)
        table = table.T

        # Add column for patients barcode
        # table['bcr_sample_barcode'] = table.index

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        return table

    def log2_transform(self, x):
        return np.log2(x + 1)


    def get_genes_list(self):
        return self.features

    def get_genes_info(self):
        if hasattr(self, "genes_info"):
            return self.genes_info

    def get_samples_list(self):
        return self.samples

    def get_network_edgelist(self):
        if hasattr(self, "network"):
            return self.network.edges()
        else:
            print(self.__class__.__str__(), "does not have network interaction data yet. (at self.network)")
            return None


class LncRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, HGNC_lncRNA_names_file_path, GENCODE_LncRNA_gtf_file_path):
        """
        :param folder_path: Path to the lncRNA expression data, downloaded from http://ibl.mdanderson.org/tanric/_design/basic/index.html
        :param lncrna_names_file_path: Path to the HGNC_RNA_long_non-coding.txt file to map ensembl gene id to lncRNA names
        """
        file_path = os.path.join(folder_path, "TCGA-rnaexpr.tsv")
        self.HGNC_lncRNA_names_path = HGNC_lncRNA_names_file_path
        self.GENCODE_LncRNA_gtf_file_path = GENCODE_LncRNA_gtf_file_path
        super().__init__(cancer_type, file_path)


    def preprocess_expression_table(self, df, columns):
        """
        Preprocess LNCRNA expression file obtained from TANRIC MDAnderson, and replace ENSEMBL gene ID to HUGO gene names (HGNC). This function overwrites the GenomicData.process_expression_table() function which processes TCGA-Assembler data.

        :param df:
        :param columns:
        :return:
        """
        lncrna_exp = df
        try:
            HGNC_lncrna_info = pd.read_table(self.HGNC_lncRNA_names_path, delimiter="\t")
            self.genes_info = HGNC_lncrna_info
        except Exception:
            raise FileNotFoundError("Needs the file RNA_long_non-coding.txt at directory external_data/HUGO_Gene_names to process lncRNA gene info")


        # Replacing ENSG Gene ID to the lncRNA gene symbol name
        lncrna_dict = self.get_lncRNA_gene_name_dict()
        # lncrna_exp['Gene_ID'] = lncrna_exp['Gene_ID'].str.replace("[.].*", "")
        lncrna_exp.replace({"Gene_ID": lncrna_dict}, inplace=True)

        # Drop NA gene rows
        lncrna_exp.dropna(axis=0, inplace=True)

        # Transpose matrix to patients rows and genes columns
        lncrna_exp.index = lncrna_exp['Gene_ID']
        lncrna_exp = lncrna_exp.T.iloc[1:, :]

        # Change index string to bcr_sample_barcode standard
        def change_patient_barcode(s):
            if "Normal" in s:
                return s[s.find('TCGA'):] + "-11A"
            elif "Tumor" in s:
                return s[s.find('TCGA'):] + "-01A"
            else:
                return s

        lncrna_exp.index = lncrna_exp.index.map(change_patient_barcode)

        return lncrna_exp

    def get_lncRNA_gene_name_dict(self):
        GENCODE_LncRNA_names = GTF.dataframe(self.GENCODE_LncRNA_gtf_file_path)
        lncrna_dict = pd.Series(GENCODE_LncRNA_names['gene_name'].values, index=GENCODE_LncRNA_names['gene_id']).to_dict()
        return lncrna_dict

    def process_starBase_miRNA_lncRNA_interactions(self, starBase_folder_path):
        self.starBase_miRNA_lncRNA_file_path = os.path.join(starBase_folder_path, "starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2018-04-26_09-10.xls")
        grn_df = pd.read_table(self.starBase_miRNA_lncRNA_file_path, header=0)

        grn_df['name'] = grn_df['name'].str.lower()
        grn_df['name'] = grn_df['name'].str.replace("-3p.*|-5p.*", "")

        self.starBase_miRNA_lncRNA_network = nx.from_pandas_dataframe(grn_df, source='name', target='geneName', create_using=nx.DiGraph())

    def get_miRNA_to_lncRNA_interactions_edgelist(self):
        return self.starBase_miRNA_lncRNA_network.edges()


class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "geneExp.txt")
        super().__init__(cancer_type, file_path)

    def process_gene_info(self, targetScan_gene_info_path, human_only=True):
        self.targetScan_gene_info_path = targetScan_gene_info_path
        self.genes_info = pd.read_table(self.targetScan_gene_info_path)

        if human_only:
            self.genes_info = self.genes_info[self.genes_info["Species ID"]==9606]

    def process_protein_coding_genes_info(self, hugo_protein_gene_names_path):
        self.hugo_protein_gene_names_path = hugo_protein_gene_names_path
        self.protein_genes_info = pd.read_table(self.hugo_protein_gene_names_path)

    def process_RegNet_gene_regulatory_network(self, grn_file_path):
        grn_df = pd.read_table(grn_file_path, header=None)
        self.regnet_grn_network = nx.from_pandas_dataframe(grn_df, source=0, target=2, create_using=nx.DiGraph())

    def get_RegNet_GRN_edgelist(self):
        return self.regnet_grn_network.edges()


class SomaticMutation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "somaticMutation_geneLevel.txt")
        super().__init__(cancer_type, file_path)


class MiRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "miRNAExp__RPM.txt")
        super().__init__(cancer_type, file_path)


    def process_target_scan(self, mirna_list, gene_symbols, targetScan_folder_path):
        self.targetScan_miR_family_info_path = os.path.join(targetScan_folder_path,"miR_Family_Info.txt")
        self.targetScan_predicted_targets_path = os.path.join(targetScan_folder_path, "Predicted_Targets_Info.default_predictions.txt")
        self.targetScan_predicted_targets_context_score_path = os.path.join(targetScan_folder_path, "Predicted_Targets_Context_Scores.default_predictions.txt")

        self.process_targetscan_mirna_family(mirna_list)
        self.process_mirna_target_interactions(mirna_list, gene_symbols)
        self.process_mirna_target_interactions_context_score(mirna_list, gene_symbols)

    def process_targetscan_mirna_family(self, mirna_list, human_only=True, incremental_group_numbering=False):
        try:
            targetScan_family_df = pd.read_table(self.targetScan_miR_family_info_path, delimiter='\t')
        except Exception:
            raise FileNotFoundError("expected TargetScan_miR_Family_Info.txt in directory mirna/TargetScan/")

        if human_only:
            targetScan_family_df = targetScan_family_df[targetScan_family_df['Species ID'] == 9606]

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19

        targetScan_family_df['MiRBase ID'] = targetScan_family_df['MiRBase ID'].str.lower()
        targetScan_family_df['MiRBase ID'] = targetScan_family_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_family_df.drop_duplicates(inplace=True)

        targetScan_family_df = targetScan_family_df[['miR family', 'MiRBase ID']]
        in_family_mirnas_list = targetScan_family_df["MiRBase ID"].tolist()
        self.mirna_family = list(targetScan_family_df["MiRBase ID"].groupby(targetScan_family_df["miR family"]))
        self.mirna_family_names = [fam[0] for fam in self.mirna_family]
        self.mirna_family = {fam[0]: fam[1].tolist() for fam in self.mirna_family}

        # Assign a unique integer number to miRNAs representing their miRNA family assignment
        self.mirna_family_assg = []
        counter = 9999
        for m in mirna_list:
            if m in in_family_mirnas_list:
                for k, v in self.mirna_family.items():
                    if m in v:
                        m_family = k
                        break
                self.mirna_family_assg.append(self.mirna_family_names.index(m_family))
            else:
                if incremental_group_numbering:
                    while counter in range(0, len(self.mirna_family_names)):
                        counter += 1
                    self.mirna_family_assg.append(counter)
                    counter += 1
                else:
                    self.mirna_family_assg.append(counter)

    def process_mirna_target_interactions(self, mirna_list, gene_symbols):
        # Load data frame from file
        try:
            targetScan_df = pd.read_table(self.targetScan_predicted_targets_path, delimiter='\t')
        except Exception:
            raise FileNotFoundError("expected TargetScan_Predicted_Targets_Info_default_predictions.txt in directory mirna/TargetScan/")

        try:
            targetScan_family_df = pd.read_table(self.targetScan_miR_family_info_path, delimiter='\t')
        except Exception:
            raise FileNotFoundError("expected TargetScan_miR_Family_Info.txt in directory mirna/TargetScan/")


        # Select only homo sapiens miRNA-target pairs
        targetScan_df = targetScan_df[targetScan_df["Species ID"] == 9606][["miR Family", "Gene Symbol"]]
        targetScan_family_df = targetScan_family_df[targetScan_family_df['Species ID'] == 9606][
            ['miR family', 'MiRBase ID']]

        # map miRBase ID names to miR Family
        targetScan_family_df.rename(columns={'miR family': 'miR Family'}, inplace=True)
        targetScan_df = pd.merge(targetScan_df, targetScan_family_df, how='inner', on="miR Family")
        targetScan_df = targetScan_df[["MiRBase ID", "Gene Symbol"]]

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        targetScan_df['MiRBase ID'] = targetScan_df['MiRBase ID'].str.lower()
        targetScan_df['MiRBase ID'] = targetScan_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_df.drop_duplicates(inplace=True)

        # Filter miRNA-target pairs to only miRNA's included in miRNA expression data, same for gene targets
        self.targetScan_df = targetScan_df[
            targetScan_df['MiRBase ID'].isin(mirna_list) & targetScan_df['Gene Symbol'].isin(gene_symbols)]

    def process_mirna_target_interactions_context_score(self, mirna_list, gene_symbols):
        # Load data frame from file
        try:
            targetScan_context_df = pd.read_table(self.targetScan_predicted_targets_context_score_path, delimiter='\t')
        except Exception:
            raise FileNotFoundError("Expected TargetScan_Predicted_Targets_Context_Scores.default_predictions.txt in directory mirna/TargetScan/")

        # Select only homo sapiens miRNA-target pairs
        targetScan_context_df = targetScan_context_df[targetScan_context_df["Gene Tax ID"] == 9606][
            ["miRNA", "Gene Symbol", "weighted context++ score percentile"]]

        # TODO Select only interactions with high context score

        # Use miRBase ID names
        targetScan_context_df.rename(columns={'miRNA': 'MiRBase ID'}, inplace=True)

        # Standardize miRNA names
        targetScan_context_df['MiRBase ID'] = targetScan_context_df['MiRBase ID'].str.lower()
        targetScan_context_df['MiRBase ID'] = targetScan_context_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_context_df.drop_duplicates(inplace=True)

        # Filter miRNA-target pairs to only miRNA's included in miRNA expression data, same for gene targets
        self.targetScan_context_df = targetScan_context_df[
            targetScan_context_df['MiRBase ID'].isin(mirna_list) & targetScan_context_df['Gene Symbol'].isin(
                gene_symbols)]

    def get_miRNA_family_group_assg(self):
        if self.mirna_family_assg is None:
            raise Exception("must first run process_target_scan(mirna_list, gene_symbols)")
        return self.mirna_family_assg

    def get_miRNA_family(self):
        if self.mirna_family is None:
            raise Exception("must first run process_target_scan(mirna_list, gene_symbols)")
        return self.mirna_family

    def get_miRNA_target_interaction(self):
        if self.targetScan_df is None:
            raise Exception("must first run process_target_scan(mirna_list, gene_symbols)")
        return self.targetScan_df

    def get_miRNA_target_interaction_context(self):
        if self.targetScan_context_df is None:
            raise Exception("must first run process_target_scan(mirna_list, gene_symbols)")
        return self.targetScan_context_df

    def get_miRNA_target_interaction_edgelist(self):
        mirna_target_interactions = self.targetScan_context_df.copy()
        mirna_target_interactions["weighted context++ score percentile"] = \
            mirna_target_interactions["weighted context++ score percentile"].apply(func=lambda x: x / 100.0)
        mirna_target_interactions.rename(columns={"weighted context++ score percentile": "weight",
                                                  "MiRBase ID": "MIR",
                                                  "Gene Symbol": "GE"}, inplace=True)
        mir_target_network = nx.from_pandas_dataframe(mirna_target_interactions,
                                                      source="MIR", target="GE", edge_attr="weight",
                                                      create_using=nx.DiGraph())
        return mir_target_network.edges(data=True)

    def get_miRNA_family_edgelist(self):
        edgelist_df = pd.DataFrame()

        for miFam in self.mirna_family.keys():
            self.mirna_family[miFam]
        #TODO finish this function


class CopyNumberVariation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "copyNumber.txt")
        super().__init__(cancer_type, file_path)


class DNAMethylation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "methylation_450.txt")
        super().__init__(cancer_type, file_path)


class ProteinExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "protein_RPPA.txt")
        super().__init__(cancer_type, file_path)

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_dataframe(HPRD_PPI, source=0, target=3,
                                          create_using=nx.DiGraph())

    def get_HPRD_PPI_network_edgelist(self):
        return self.HPRD_PPI_network.edges()

    def process_STRING_PPI_network(self, ppi_data_file_path):
        pass


