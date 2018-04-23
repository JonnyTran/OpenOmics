import os

from definitions import ROOT_DIR

import numpy as np
import pandas as pd


class GenomicData:
    def __init__(self, cancer_type, file_path, columns="GeneSymbol|TCGA", log2_transform=True):
        self.cancer_type = cancer_type

        self.data = self.preprocess_expression_table(pd.read_table(file_path, sep="\t"), columns)

        if log2_transform:
            self.data = self.data.applymap(self.log2_transform)

        # Save samples and features for this omics data
        self.samples = self.data.index
        self.features = self.data.columns.tolist()
        # self.features.remove("bcr_sample_barcode")

    def preprocess_expression_table(self, df, columns):
        """
        Download
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


    def get_samples_list(self):
        return self.samples


class LncRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "TCGA-rnaexpr.tsv")
        super().__init__(cancer_type, file_path)

    def preprocess_expression_table(self, df, columns):
        lncrna_exp = df

        lncrna_names = pd.read_table(
            os.path.join(ROOT_DIR, "data/tcga-assembler/HGNC_RNA_long_non-coding.txt"),
            delimiter="\t")
        lncrna_dict = pd.Series(lncrna_names.symbol.values, index=lncrna_names.ensembl_gene_id).to_dict()

        # Replacing ENSG Gene ID to the lncRNA symbol name
        lncrna_exp['Gene_ID'] = lncrna_exp['Gene_ID'].str.replace("[.].*", "")
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


class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "geneExp.txt")
        super().__init__(cancer_type, file_path)


class SNP(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "somaticMutation_geneLevel.txt")
        super().__init__(cancer_type, file_path)


class miRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "miRNAExp__RPM.txt")
        super().__init__(cancer_type, file_path)

    def process_target_scan(self, mirna_list, gene_symbols):
        self.process_targetscan_mirna_family(mirna_list)
        self.process_mirna_target_interactions(mirna_list, gene_symbols)
        self.process_mirna_target_interactions_context_score(mirna_list, gene_symbols)

    def process_targetscan_mirna_family(self, mirna_list, incremental_group_numbering=False):
        targetScan_family_df = pd.read_table(
            os.path.join(ROOT_DIR, 'data/external/TargetScan/TargetScan_miR_Family_Info.txt'),
            delimiter='\t')
        targetScan_family_df = targetScan_family_df[targetScan_family_df['Species ID'] == 9606]
        targetScan_family_df['MiRBase ID'] = targetScan_family_df['MiRBase ID'].str.lower()
        targetScan_family_df['MiRBase ID'] = targetScan_family_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_family_df.drop_duplicates(inplace=True)
        targetScan_family_df = targetScan_family_df[['miR family', 'MiRBase ID']]
        in_family_mirnas_list = targetScan_family_df["MiRBase ID"].tolist()
        self.mirna_family = list(targetScan_family_df["MiRBase ID"].groupby(targetScan_family_df["miR family"]))
        self.mirna_family_names = [fam[0] for fam in self.mirna_family]
        self.mirna_family = {fam[0]: fam[1].tolist() for fam in self.mirna_family}

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
        targetScan_df = pd.read_table(
            os.path.join(ROOT_DIR,
                         'data/external/TargetScan/TargetScan_Predicted_Targets_Info_default_predictions.tsv'),
            delimiter='\t')
        targetScan_family_df = pd.read_table(
            os.path.join(ROOT_DIR, 'data/external/TargetScan/TargetScan_miR_Family_Info.txt'),
            delimiter='\t')

        # Select only homo sapiens miRNA-target pairs
        targetScan_df = targetScan_df[targetScan_df["Species ID"] == 9606][["miR Family", "Gene Symbol"]]
        targetScan_family_df = targetScan_family_df[targetScan_family_df['Species ID'] == 9606][
            ['miR family', 'MiRBase ID']]

        # Use miRBase ID names
        targetScan_family_df.rename(columns={'miR family': 'miR Family'}, inplace=True)
        targetScan_df = pd.merge(targetScan_df, targetScan_family_df, how='inner', on="miR Family")
        targetScan_df = targetScan_df[["MiRBase ID", "Gene Symbol"]]

        # Standardize miRNA names
        targetScan_df['MiRBase ID'] = targetScan_df['MiRBase ID'].str.lower()
        targetScan_df['MiRBase ID'] = targetScan_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_df.drop_duplicates(inplace=True)

        # Filter miRNA-target pairs to only miRNA's included in miRNA expression data, same for gene targets
        self.targetScan_df = targetScan_df[
            targetScan_df['MiRBase ID'].isin(mirna_list) & targetScan_df['Gene Symbol'].isin(gene_symbols)]

    def process_mirna_target_interactions_context_score(self, mirna_list, gene_symbols):
        # Load data frame from file
        targetScan_context_df = pd.read_table(
            os.path.join(ROOT_DIR,
                         'data/external/TargetScan/TargetScan_Predicted_Targets_Context_Scores.default_predictions.txt'),
            delimiter='\t')

        # Select only homo sapiens miRNA-target pairs
        targetScan_context_df = targetScan_context_df[targetScan_context_df["Gene Tax ID"] == 9606][
            ["miRNA", "Gene Symbol"]]

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
        return self.mirna_family_assg

    def get_miRNA_target_interaction(self):
        return self.targetScan_df

    def get_miRNA_target_interaction_context(self):
        return self.targetScan_context_df



class CopyNumberVariation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "copyNumber.txt")
        super().__init__(cancer_type, file_path)


class DNAMethylation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "")
        super().__init__(cancer_type, file_path)


class ProteinExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "protein_RPPA.txt")
        super().__init__(cancer_type, file_path)


if __name__ == '__main__':
    folder_path = "/data/tcga-assembler/LUAD/lncrna/"
    lncRNA_expression = LncRNAExpression(cancer_type="LUAD", folder_path=ROOT_DIR + folder_path)
