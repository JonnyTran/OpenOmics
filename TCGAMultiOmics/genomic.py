import os

import networkx as nx
import numpy as np
import pandas as pd

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

        TANRIC LNCRNA expression values are log2 transformed

        :param df:
        :param columns:
        :return:
        """
        lncrna_exp = df
        try:
            HGNC_lncrna_info = pd.read_table(self.HGNC_lncRNA_names_path, delimiter="\t", usecols=['symbol', 'ensembl_gene_id', 'name', 'location'])
            self.HGNC_lncrna_info = HGNC_lncrna_info
            self.HGNC_lncrna_info.index = self.HGNC_lncrna_info["symbol"]
        except Exception:
            raise FileNotFoundError("Needs the file RNA_long_non-coding.txt at directory external_data/HUGO_Gene_names to process lncRNA gene info")


        # Replacing ENSG Gene ID to the lncRNA gene symbol name
        ensembl_id_to_gene_name, ensembl_id_to_transcript_id = self.get_GENCODE_lncRNA_gene_name_dict()
        hgnc_lncrna_dict = self.get_HUGO_lncRNA_gene_name_dict()
        lncrna_exp['Gene_ID'] = lncrna_exp['Gene_ID'].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end

        # Preprocess genes info
        self.preprocess_genes_info(lncrna_exp['Gene_ID'],     ensembl_id_to_gene_name, ensembl_id_to_transcript_id, hgnc_lncrna_dict)

        lncrna_exp.replace({"Gene_ID": ensembl_id_to_gene_name}, inplace=True)
        lncrna_exp.replace({"Gene_ID": hgnc_lncrna_dict}, inplace=True)

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


    def get_GENCODE_lncRNA_gene_name_dict(self):
        GENCODE_LncRNA_names = GTF.dataframe(self.GENCODE_LncRNA_gtf_file_path)

        GENCODE_LncRNA_names['gene_id'] = GENCODE_LncRNA_names['gene_id'].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end

        ensembl_id_to_gene_name = pd.Series(GENCODE_LncRNA_names['gene_name'].values, index=GENCODE_LncRNA_names['gene_id']).to_dict()

        ensembl_id_to_transcript_id = pd.Series(GENCODE_LncRNA_names['transcript_id'].values, index=GENCODE_LncRNA_names['gene_id']).to_dict()

        return ensembl_id_to_gene_name, ensembl_id_to_transcript_id

    def get_HUGO_lncRNA_gene_name_dict(self):
        lncrna_dict = pd.Series(self.HGNC_lncrna_info['symbol'].values,
                                index=self.HGNC_lncrna_info['ensembl_gene_id']).to_dict()
        return lncrna_dict

    def process_starBase_miRNA_lncRNA_interactions(self, starBase_folder_path):
        self.starBase_miRNA_lncRNA_file_path = os.path.join(starBase_folder_path, "starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2018-04-26_09-10.xls")
        grn_df = pd.read_table(self.starBase_miRNA_lncRNA_file_path, header=0)

        grn_df['name'] = grn_df['name'].str.lower()
        grn_df['name'] = grn_df['name'].str.replace("-3p.*|-5p.*", "")

        self.starBase_miRNA_lncRNA_network = nx.from_pandas_dataframe(grn_df, source='geneName', target='name', create_using=nx.DiGraph())

    def get_starBase_lncRNA_miRNA_interactions_edgelist(self):
        return self.starBase_miRNA_lncRNA_network.edges()

    def process_lncRNome_miRNA_binding_sites(self, lncRNome_folder_path):
        self.lnRNome_miRNA_binding_sites_path = os.path.join(lncRNome_folder_path, "miRNA_binding_sites.txt")

        df = pd.read_table(self.lnRNome_miRNA_binding_sites_path, header=0)

        df['Binding miRNAs'] = df['Binding miRNAs'].str.lower()
        df['Binding miRNAs'] = df['Binding miRNAs'].str.replace("-3p.*|-5p.*", "")

        self.lncRNome_miRNA_binding_sites_network = nx.from_pandas_dataframe(df, source='Gene Name', target='Binding miRNAs', create_using=nx.DiGraph())

    def get_lncRNome_miRNA_binding_sites_edgelist(self):
        return self.lncRNome_miRNA_binding_sites_network.edges()

    def preprocess_genes_info(self, genes_list, ensembl_id_to_gene_name, ensembl_id_to_transcript_id, hugo_lncrna_dict):
        self.gene_info = pd.DataFrame(index=genes_list)
        self.gene_info.index.name = "genes list"

        self.gene_info["Gene Name"] = self.gene_info.index.map(ensembl_id_to_gene_name)
        self.gene_info["Gene Name"].fillna({"ensembl id": hugo_lncrna_dict}, inplace=True)

        self.gene_info["Transcript id"] = self.gene_info.index.map(ensembl_id_to_transcript_id)

        self.gene_info.index = genes_list


    def process_lncRNome_gene_info(self, lncRNome_folder_path):
        self.lnRNome_genes_info_path = os.path.join(lncRNome_folder_path, "general_information.txt")

        self.lnRNome_genes_info = pd.read_table(self.lnRNome_genes_info_path, header=0, usecols=["Gene Name", "Transcript Name", "Transcript Type", "Location", "Strand"])


    def process_NONCODE_func_annotation(self, noncode_folder_path):
        self.noncode_source_path = os.path.join(noncode_folder_path, "NONCODEv5_source")
        self.noncode_transcript2gene_path = os.path.join(noncode_folder_path, "NONCODEv5_Transcript2Gene")
        self.noncode_func_path = os.path.join(noncode_folder_path, "NONCODEv5_human.func")

        # Import tables
        source_df = pd.read_table(self.noncode_source_path, header=None)
        source_df.columns = ["NONCODE Transcript ID", "name type", "Gene ID"]

        transcript2gene_df = pd.read_table(self.noncode_transcript2gene_path, header=None)
        transcript2gene_df.columns = ["NONCODE Transcript ID", "NONCODE Gene ID"]

        self.noncode_func_df = pd.read_table(self.noncode_func_path, header=None)
        self.noncode_func_df.columns = ["NONCODE Gene ID", "GO terms"]
        self.noncode_func_df.index = self.noncode_func_df["NONCODE Gene ID"]

        # Convert to NONCODE transcript ID for the functional annotattion data
        self.noncode_func_df["NONCODE Transcript ID"] = self.noncode_func_df.index.map(pd.Series(transcript2gene_df['NONCODE Transcript ID'].values, index=transcript2gene_df['NONCODE Gene ID']).to_dict())

        # Convert NONCODE transcript ID to gene names
        source_df = source_df[source_df["name type"] == "NAME"]

        self.noncode_func_df["Gene Name"] = self.noncode_func_df["NONCODE Transcript ID"].map(pd.Series(source_df['Gene ID'].values,index=source_df['NONCODE Transcript ID']).to_dict())


    def get_genes_info(self):
        if ~hasattr(self, "genes_info_processed") or self.genes_info_processed == False:
            # self.gene_info = pd.merge(self.gene_info, self.HGNC_lncrna_info.groupby("symbol").first(), how="left", left_on="Gene Name", right_on="symbol")
            self.gene_info.index.name = "symbol"
            self.gene_info = self.gene_info.join(self.HGNC_lncrna_info.groupby("symbol").first(), on="symbol",
                                                 how="left")

            # self.gene_info = pd.merge(self.gene_info, self.lnRNome_genes_info.groupby("Gene Name").first(), how="left", left_on="Gene Name", right_on="Gene Name")
            self.gene_info.index.name = "Gene Name"
            self.gene_info = self.gene_info.join(self.lnRNome_genes_info.groupby("Gene Name").first(), on="Gene Name",
                                                 how="left")

            # self.gene_info = pd.merge(self.gene_info, self.noncode_func_df.groupby("Gene Name").first(), how="left", left_on="Gene Name", right_on="Gene Name")
            self.gene_info.index = self.gene_info["Gene Name"]
            self.gene_info = self.gene_info.join(self.noncode_func_df.groupby("Gene Name").first(), on="Gene Name",
                                                 how="left")

            self.gene_info.index = self.get_genes_list()
            self.genes_info_processed = True
        return self.gene_info



class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "geneExp.txt")
        super().__init__(cancer_type, file_path)

    def process_targetScan_gene_info(self, targetScan_gene_info_path, human_only=True):
        self.targetScan_gene_info_path = targetScan_gene_info_path
        self.targetScan_genes_info = pd.read_table(self.targetScan_gene_info_path, usecols=["Transcript ID", "Gene ID", "Species ID", "Gene symbol", "Gene description", "3P-seq tags"])

        self.targetScan_genes_info["Gene description"] = self.targetScan_genes_info["Gene description"].str.replace(" \[.*\]","")

        if human_only:
            self.targetScan_genes_info = self.targetScan_genes_info[self.targetScan_genes_info["Species ID"] == 9606]
        self.targetScan_genes_info.drop(columns=["Species ID"], inplace=True)


    def process_protein_coding_genes_info(self, hugo_protein_gene_names_path):
        self.hugo_protein_gene_names_path = hugo_protein_gene_names_path
        self.hugo_protein_genes_info = pd.read_table(self.hugo_protein_gene_names_path, usecols=["symbol", "locus_type", "gene_family", "gene_family_id", "location"])


    def process_RegNet_gene_regulatory_network(self, grn_file_path):
        grn_df = pd.read_table(grn_file_path, header=None)

        # Since RegNet GRN contains miRNA and TF regulatory interactions
        # hsa-miR-* microRNA gene names will be mapped to hsa-mir-*
        grn_df[0] = grn_df[0].map(lambda x: x.lower() if ("hsa-miR" in x) else x)

        self.regnet_grn_network = nx.from_pandas_dataframe(grn_df, source=0, target=2, create_using=nx.DiGraph())


    def get_RegNet_GRN_edgelist(self):
        return self.regnet_grn_network.edges()

    def get_genes_info(self):
        gene_info = pd.DataFrame(index=self.get_genes_list())

        gene_info.index.name = "Gene symbol"
        gene_info = gene_info.join(self.targetScan_genes_info.groupby("Gene symbol").first(), on="Gene symbol", how="left")

        gene_info.index.name = "symbol"
        gene_info = gene_info.join(self.hugo_protein_genes_info.groupby("symbol").first(), on="symbol", how="left")
        return gene_info


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

        targetScan_family_df = targetScan_family_df[['miR family', 'MiRBase ID', 'Seed+m8', 'Mature sequence', 'Family Conservation?', 'MiRBase Accession']]

        targetScan_family_df['MiRBase ID'] = targetScan_family_df['MiRBase ID'].astype(str)

        self.targetScan_family_df = targetScan_family_df
        self.targetScan_family_df.index = self.targetScan_family_df['MiRBase ID']


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
        targetScan_family_df = targetScan_family_df[targetScan_family_df['Species ID'] == 9606][['miR family', 'MiRBase ID']]

        # map miRBase ID names to miR Family
        targetScan_family_df.rename(columns={'miR family': 'miR Family'}, inplace=True)
        targetScan_df = pd.merge(targetScan_df, targetScan_family_df, how='inner', on="miR Family")
        targetScan_df = targetScan_df[["MiRBase ID", "Gene Symbol"]]

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        targetScan_df['MiRBase ID'] = targetScan_df['MiRBase ID'].str.lower()
        targetScan_df['MiRBase ID'] = targetScan_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_df.drop_duplicates(inplace=True)


    def process_mirna_target_interactions_context_score(self, mirna_list, gene_symbols):
        # Load data frame from file
        try:
            targetScan_context_df = pd.read_table(self.targetScan_predicted_targets_context_score_path, delimiter='\t')
        except Exception:
            raise FileNotFoundError("Expected TargetScan_Predicted_Targets_Context_Scores.default_predictions.txt in directory mirna/TargetScan/")

        # Select only homo sapiens miRNA-target pairs
        targetScan_context_df = targetScan_context_df[targetScan_context_df["Gene Tax ID"] == 9606][
            ["miRNA", "Gene Symbol", "weighted context++ score percentile"]]

        # Use miRBase ID names
        targetScan_context_df.rename(columns={'miRNA': 'MiRBase ID'}, inplace=True)

        # Standardize miRNA names
        targetScan_context_df['MiRBase ID'] = targetScan_context_df['MiRBase ID'].str.lower()
        targetScan_context_df['MiRBase ID'] = targetScan_context_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        targetScan_context_df.drop_duplicates(inplace=True)


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


    def get_genes_info(self):
        gene_info = pd.DataFrame(index=self.get_genes_list())

        gene_info.index.name = "MiRBase ID"
        gene_info = gene_info.join(self.targetScan_family_df.groupby("MiRBase ID").first(), on="MiRBase ID",how="left")

        return gene_info


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


