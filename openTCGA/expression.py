import os
from collections import OrderedDict

import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.UniProt import GOA
from pandas import Series

from openTCGA.utils import GTF


class ExpressionData:
    def __init__(self, cohort_name, file_path, columns, key,
                 import_sequences="longest", replace_U2T=True,
                 transposed_table=True, log2_transform=False):
        """
        .. class:: ExpressionData
        An abstract class that handles importing of expression data tables while providing indices to the TCGA
        samples and gene name to the expressions.
            Args:
                cohort_name (str): the cohort code name string
                file_path (str): Path of the table file to import
                columns (str): column names to import from the table. Columns names imported are string match, separated by "|"
                transposed_table (list): ["longest", "shortest", "multi"], If True, perform preprocessing steps for the table data obtained from TCGA-Assembler tool. If False, import a pandas table as-is with bcr_sample_barcode for row index, and gene names as columns
                log2_transform (bool): Whether to log2 transform the expression values
        """
        self.cohort_name = cohort_name
        self.import_sequences = import_sequences
        self.replace_U2T = replace_U2T

        if os.path.isfile(file_path) and os.path.exists(file_path):
            table = pd.read_table(file_path)
        else:
            raise FileNotFoundError(file_path)

        if transposed_table:
            self.expression = self.preprocess_table(table, columns, key)

        if log2_transform:
            self.expression = self.expression.applymap(self.log2_transform)

        # Save samples and features for this omics data
        self.samples = self.expression.index
        self.features = self.expression.columns.tolist()
        # self.features.remove("bcr_sample_barcode")

    def preprocess_table(self, table:pd.DataFrame, columns:str, key:str):
        """
        This function preprocesses the expression table files where columns are samples and rows are gene/transcripts
        :param df:
        :param columns:
        :return:
        """
        # Filter columns
        table = table.filter(regex=columns)

        # Cut TCGA column names to sample barcode, discarding aliquot info
        table = table.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x, inplace=False)

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        # Drop NA geneID rows
        table.dropna(axis=0, inplace=True)

        # Remove entries with unknown geneID
        table = table[table[key] != '?']

        # Transpose dataframe to patient rows and geneID columns
        table.index = table[key]
        table.drop([key], axis=1, inplace=True)
        table = table.T

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        return table

    def log2_transform(self, x):
        return np.log2(x + 1)

    def drop_genes(self, genes_to_drop):
        self.expression.drop(genes_to_drop, axis=1, inplace=True)
        for gene in genes_to_drop:
            self.features.remove(gene)

    def get_genes_list(self):
        return self.features

    def get_genes_info(self):
        if hasattr(self, "genes_info"):
            return self.genes_info

    def get_samples_list(self):
        return self.samples

    def get_network_edgelist(self):
        if hasattr(self, "network"):
            return self.network.edges(data=True)
        else:
            print(self.__class__.__str__(), "does not have network interaction data yet. (at self.network)")
            return None


class LncRNAExpression(ExpressionData):
    def __init__(self, cohort_name, file_path, columns="Gene_ID|TCGA", key="Gene_ID",
                 HGNC_lncRNA_names_file_path=None, GENCODE_folder_path=None, external_data_path=None):
        """
        :param file_path: Path to the lncRNA expression data, downloaded from http://ibl.mdanderson.org/tanric/_design/basic/index.html
        """
        self.HGNC_lncRNA_names_path = HGNC_lncRNA_names_file_path
        self.GENCODE_LncRNA_gtf_file_path = os.path.join(GENCODE_folder_path, "gencode.v29.long_noncoding_RNAs.gtf")
        self.GENCODE_LncRNA_sequence_file_path = os.path.join(GENCODE_folder_path, "gencode.v29.lncRNA_transcripts.fa")
        self.external_data_path = external_data_path
        super().__init__(cohort_name, file_path, columns=columns, key=key)

    def preprocess_table(self, df, columns, key):
        """
        Preprocess LNCRNA expression file obtained from TANRIC MDAnderson, and replace ENSEMBL gene ID to HUGO gene names (HGNC). This function overwrites the GenomicData.process_expression_table() function which processes TCGA-Assembler data.

        TANRIC LNCRNA expression values are log2 transformed
        """
        lncrna_exp = df

        # Replacing ENSG Gene ID to the lncRNA gene symbol name
        lncrna_exp[key] = lncrna_exp[key].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
        lncrna_exp = lncrna_exp[~lncrna_exp[key].duplicated(keep='first')] # Remove duplicate genes

        # Preprocess genes info
        gencode_LncRNA_info, ensembl_gene_id_to_gene_name = self.get_GENCODE_lncRNA_gene_name_dict()
        lncipedia_lncrna_dict = self.get_lncipedia_gene_id_to_name_dict()
        lncBase_gene_id_to_name_dict = self.get_lncBase_gene_id_to_name_dict()

        hgnc_lncrna_dict = self.get_HUGO_lncRNA_gene_name_dict()
        ensembl_gene_ids = lncrna_exp[key]
        self.preprocess_genes_info(ensembl_gene_ids, gencode_LncRNA_info,
                                   ensembl_gene_id_to_gene_name,
                                   hgnc_lncrna_dict)

        # Convert ensembl gene IDs to known gene names
        print("Unmatched lncRNAs", lncrna_exp[key].str.startswith("ENSG").sum())

        lncrna_exp.replace({key: ensembl_gene_id_to_gene_name}, inplace=True)
        print("Unmatched lncRNAs after gencode:", lncrna_exp['Gene_ID'].str.startswith("ENSG").sum())

        lncrna_exp.replace({key: lncBase_gene_id_to_name_dict}, inplace=True)
        print("Unmatched lncRNAs after lncBase:", lncrna_exp['Gene_ID'].str.startswith("ENSG").sum())

        lncrna_exp.replace({key: hgnc_lncrna_dict}, inplace=True)
        print("Unmatched lncRNAs after HGNC:", lncrna_exp['Gene_ID'].str.startswith("ENSG").sum())

        lncrna_exp.replace({key: lncipedia_lncrna_dict}, inplace=True)
        print("Unmatched lncRNAs after lncipedia:", lncrna_exp['Gene_ID'].str.startswith("ENSG").sum())

        # Drop NA gene rows
        lncrna_exp.dropna(axis=0, inplace=True)

        # Transpose matrix to patients rows and genes columns
        lncrna_exp.index = lncrna_exp[key]
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

    def get_lncBase_gene_id_to_name_dict(self):
        table = pd.read_table(os.path.join(self.external_data_path, "lncBase/LncBasev2_download.csv"))
        lncBase_gene_id_to_name_dict = pd.Series(table["geneName"].values,
                                          index=table["geneId"]).to_dict()
        return lncBase_gene_id_to_name_dict

    def get_lncipedia_gene_id_to_name_dict(self):
        lncipedia_names = GTF.dataframe(
            os.path.join(self.external_data_path, "LNCipedia/lncipedia_5_0_hg19 (copy).gtf"))
        lncipedia_names = lncipedia_names[
            lncipedia_names["gene_alias_1"].notnull() & lncipedia_names["gene_alias_2"].notnull()]
        lncipedia_names = lncipedia_names[lncipedia_names["gene_alias_1"].str.startswith("ENSG")]
        lncipedia_names = lncipedia_names[~lncipedia_names["gene_alias_2"].str.startswith("ENSG")]
        lncipedia_lncrna_dict = pd.Series(lncipedia_names["gene_alias_2"].values,
                                          index=lncipedia_names["gene_alias_1"]).to_dict()
        return lncipedia_lncrna_dict

    def preprocess_genes_info(self, ensembl_gene_id, GENCODE_LncRNA_info, ensembl_id_to_gene_name, hgnc_lncrna_dict):
        self.gene_info = pd.DataFrame(index=ensembl_gene_id)
        self.gene_info.index.name = "ensembl_gene_id"

        self.gene_info["Gene ID"] = self.gene_info.index
        self.gene_info["Gene Name"] = self.gene_info.index.map(ensembl_id_to_gene_name)
        self.gene_info["HGNC Gene Name"] = self.gene_info.index.map(hgnc_lncrna_dict)


        self.gene_info["Transcript id"] = self.gene_info.index.map(GENCODE_LncRNA_info[
            GENCODE_LncRNA_info["transcript_id"].notnull()].groupby('gene_id')["transcript_id"].apply(
            lambda x: "|".join(x.unique())).to_dict())
        self.gene_info["Transcript name"] = self.gene_info.index.map(
            GENCODE_LncRNA_info[
                GENCODE_LncRNA_info["transcript_name"].notnull()].groupby('gene_id')["transcript_name"].apply(
                lambda x: "|".join(x.unique())).to_dict())

        self.gene_info["Transcript type"] = self.gene_info.index.map(
            GENCODE_LncRNA_info[GENCODE_LncRNA_info["transcript_type"].notnull()].groupby('gene_id')["transcript_type"].apply(
                lambda x: "|".join(x.unique())).to_dict())

        self.gene_info["tag"] = self.gene_info.index.map(
            GENCODE_LncRNA_info[GENCODE_LncRNA_info["tag"].notnull()].groupby('gene_id')["tag"].apply(
                lambda x: "|".join(x.unique())).to_dict())

        self.gene_info["Chromosome"] = self.gene_info.index.map(pd.Series(GENCODE_LncRNA_info['seqname'].values,
                                                                               index=GENCODE_LncRNA_info['gene_id']).to_dict())
        self.gene_info["start"] = self.gene_info.index.map(pd.Series(GENCODE_LncRNA_info['start'].values,
                                                                          index=GENCODE_LncRNA_info[
                                                                              'gene_id']).to_dict()).astype(np.float64)
        self.gene_info["end"] = self.gene_info.index.map(pd.Series(GENCODE_LncRNA_info['end'].values,
                                                                     index=GENCODE_LncRNA_info[
                                                                         'gene_id']).to_dict()).astype(np.float64)
        self.gene_info["Strand"] = self.gene_info.index.map(pd.Series(GENCODE_LncRNA_info['strand'].values,
                                                                   index=GENCODE_LncRNA_info[
                                                                       'gene_id']).to_dict())

        self.gene_info["locus_type"] = self.gene_info.index.map(pd.Series(GENCODE_LncRNA_info['gene_type'].values,
                                                                     index=GENCODE_LncRNA_info[
                                                                         'gene_id']).to_dict())


        # Merge GENCODE transcript sequence data
        self.gene_info["Transcript sequence"] = self.gene_info["Gene Name"].map(
            self.get_GENCODE_lncRNA_sequence_data(self.import_sequences, self.replace_U2T))


    def get_GENCODE_lncRNA_gene_name_dict(self):
        GENCODE_LncRNA_info = GTF.dataframe(self.GENCODE_LncRNA_gtf_file_path)
        # print(GENCODE_LncRNA_info.columns)
        GENCODE_LncRNA_info['gene_id'] = GENCODE_LncRNA_info['gene_id'].str.replace("[.].*", "")  # TODO Removing .# ENGS gene version number at the end
        GENCODE_LncRNA_info['transcript_id'] = GENCODE_LncRNA_info['transcript_id'].str.replace("[.].*", "")

        ensembl_id_to_gene_name = pd.Series(GENCODE_LncRNA_info['gene_name'].values,
                                            index=GENCODE_LncRNA_info['gene_id']).to_dict()

        return GENCODE_LncRNA_info, ensembl_id_to_gene_name

    def get_HUGO_lncRNA_gene_name_dict(self):
        try:
            HGNC_lncrna_info = pd.read_table(self.HGNC_lncRNA_names_path, delimiter="\t",
                                             usecols=['symbol', 'locus_type', 'ensembl_gene_id', 'name', 'location'])
            self.HGNC_lncrna_info = HGNC_lncrna_info
            self.HGNC_lncrna_info.index = self.HGNC_lncrna_info["ensembl_gene_id"]
        except Exception:
            raise FileNotFoundError("Needs the file RNA_long_non-coding.txt at directory external_data/HUGO_Gene_names to process lncRNA gene info")

        lncrna_dict = pd.Series(self.HGNC_lncrna_info['symbol'].values,
                                index=self.HGNC_lncrna_info['ensembl_gene_id']).to_dict()
        return lncrna_dict

    def process_starBase_miRNA_lncRNA_interactions(self, starBase_folder_path):
        self.starBase_miRNA_lncRNA_file_path = os.path.join(starBase_folder_path,
                                                            "starBase_Human_Pan-Cancer_miRNA-LncRNA_Interactions2018-04-26_09-10.xls")

    def process_starBase_lncRNA_RNA_interactions(self, starBase_folder_path):
        self.starBase_lncRNA_RNA_interactions_file_path = os.path.join(starBase_folder_path,
                                                            "starbase_3.0_lncrna_rna_interactions.csv")

    def process_LncReg_lncRNA_RNA_regulatory_interactions(self, LncReg_folder_path):
        self.LncReg_RNA_regulatory_file_path = os.path.join(LncReg_folder_path, "data.xlsx")

    def process_lncBase_miRNA_lncRNA_interactions(self, lncBase_folder_path):
        self.lncBase_interactions_file_path = os.path.join(lncBase_folder_path, "LncBasev2_download.csv")
        self.lncBase_predicted_interactions_file_path = os.path.join(lncBase_folder_path, "lncBaseV2_predicted_human_data.csv")

    def process_NPInter_ncRNA_RNA_regulatory_interactions(self, NPInter_folder_path):
        self.NPInter_interactions_file_path = os.path.join(NPInter_folder_path, "interaction_NPInter[v3.0].txt")
        self.NPInter_interactions_old_file_path = os.path.join(NPInter_folder_path, "interaction_NPInter[v2.0].txt")

    def process_lncRNome_miRNA_binding_sites(self, lncRNome_folder_path):
        self.lnRNome_miRNA_binding_sites_path = os.path.join(lncRNome_folder_path, "miRNA_binding_sites.txt")

    def get_starBase_miRNA_lncRNA_interactions_edgelist(self, data=True, rename_dict=None):
        grn_df = pd.read_table(self.starBase_miRNA_lncRNA_file_path, header=0)

        grn_df['name'] = grn_df['name'].str.lower()
        grn_df['name'] = grn_df['name'].str.replace("-3p.*|-5p.*", "")

        self.starBase_miRNA_lncRNA_network = nx.from_pandas_edgelist(grn_df, source='name', target='geneName',
                                                                     create_using=nx.DiGraph())
        if rename_dict is not None:
            self.starBase_miRNA_lncRNA_network = nx.relabel_nodes(self.starBase_miRNA_lncRNA_network, rename_dict)
        return self.starBase_miRNA_lncRNA_network.edges(data=data)

    def get_starBase_lncRNA_RNA_interactions(self, data=True, rename_dict=None):
        df = pd.read_table(self.starBase_lncRNA_RNA_interactions_file_path, header=0, sep=",",
              usecols=["geneID", "geneName", "geneType", "pairGeneID", "pairGeneName",
                      "pairGeneType", "interactionNum", 'expNum', "FreeEnergy"])

        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[df["pairGeneType"] == "miRNA"][
            "pairGeneName"].str.lower()
        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[df["pairGeneType"] == "miRNA"][
            "pairGeneName"].str.replace("-3p.*|-5p.*", "")


        self.starBase_lncRNA_RNA_network = nx.from_pandas_edgelist(df, source='geneName', target='pairGeneName',
                                                                   edge_attr=["interactionNum"],
                                                                     create_using=nx.DiGraph())
        if rename_dict is not None:
            self.starBase_lncRNA_RNA_network = nx.relabel_nodes(self.starBase_lncRNA_RNA_network, rename_dict)
        return self.starBase_lncRNA_RNA_network.edges(data=data)

    def get_LncReg_lncRNA_RNA_regulatory_interactions(self, data=True):
        table = pd.read_excel(self.LncReg_RNA_regulatory_file_path)

        table = table[table["species"] == "Homo sapiens"]
        table.loc[table["B_category"] == "miRNA", "B_name_in_paper"] = table[table["B_category"] == "miRNA"][
            "B_name_in_paper"].str.replace("-3p.*|-5p.*", "")
        table.loc[table["B_category"] == "miRNA", "B_name_in_paper"] = table[table["B_category"] == "miRNA"][
            "B_name_in_paper"].str.replace("MIR", "hsa-mir-")
        table.loc[table["B_category"] == "miRNA", "B_name_in_paper"] = table[table["B_category"] == "miRNA"][
            "B_name_in_paper"].str.replace("let-", "hsa-let-")

        LncReg_lncRNA_RNA_network = nx.from_pandas_edgelist(table, source='A_name_in_paper', target='B_name_in_paper',
                                                               edge_attr=["relationship", "mechanism", "pmid"],
                                                               create_using=nx.DiGraph())
        return LncReg_lncRNA_RNA_network.edges(data=data)


    def get_lncBase_miRNA_lncRNA_interactions_edgelist(self, tissue=None, data=True, rename_dict=None):
        lncbase_df = pd.read_table(self.lncBase_interactions_file_path)

        lncbase_df = lncbase_df[lncbase_df["species"] == "Homo sapiens"]
        lncbase_df["mirna"] = lncbase_df["mirna"].str.lower()
        lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("-3p.*|-5p.*", "")

        if tissue is not None:
            lncbase_df = lncbase_df[lncbase_df["tissue"] == tissue]

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source='mirna', target='geneId',
                                                               edge_attr=["tissue", "positive_negative"],
                                                               create_using=nx.DiGraph())
        if rename_dict is not None:
            lncBase_lncRNA_miRNA_network = nx.relabel_nodes(lncBase_lncRNA_miRNA_network, rename_dict)

        return lncBase_lncRNA_miRNA_network.edges(data=data)

    def get_lncBase_miRNA_lncRNA_predicted_interactions_edgelist(self, data=True, rename_dict=None):
        records = []
        for record in SeqIO.parse(
                self.lncBase_predicted_interactions_file_path,
                "fasta"):
            records.append(record.id.split(","))
        lncbase_df = pd.DataFrame(records, columns=["Transcript_ID", "geneId", "mirna", "miTG-score"])

        lncbase_df["miTG-score"] = lncbase_df["miTG-score"].astype(float)
        lncbase_df = lncbase_df[lncbase_df["miTG-score"] >= 0.9]

        lncbase_df["geneId"] = lncbase_df["geneId"].str.replace("(\(.+|\))", "")
        lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("(\(.+|\))", "")
        lncbase_df["mirna"] = lncbase_df["mirna"].str.lower()
        lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("-3p.*|-5p.*", "")

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source='mirna', target='geneId',
                                                               edge_attr=["miTG-score"],
                                                               create_using=nx.DiGraph())
        if rename_dict is not None:
            lncBase_lncRNA_miRNA_network = nx.relabel_nodes(lncBase_lncRNA_miRNA_network, rename_dict)

        return lncBase_lncRNA_miRNA_network.edges(data=data)

    def get_lncRNome_miRNA_binding_sites_edgelist(self, data=True, rename_dict=None):
        df = pd.read_table(self.lnRNome_miRNA_binding_sites_path, header=0)

        df['Binding miRNAs'] = df['Binding miRNAs'].str.lower()
        df['Binding miRNAs'] = df['Binding miRNAs'].str.replace("-3p.*|-5p.*", "")

        lncRNome_miRNA_binding_sites_network = nx.from_pandas_edgelist(df, source='Gene Name',
                                                                            target='Binding miRNAs',
                                                                            edge_attr=["miRNA Interaction Site",
                                                                                       "Transcript ID"],
                                                                            create_using=nx.DiGraph())
        if rename_dict is not None:
            lncRNome_miRNA_binding_sites_network = nx.relabel_nodes(lncRNome_miRNA_binding_sites_network, rename_dict)

        return lncRNome_miRNA_binding_sites_network.edges(data=data)


    def get_NPInter_ncRNA_RNA_regulatory_interaction_edgelist(self, use_latest=True, data=True, rename_dict=None):
        if use_latest:
            file_path = self.NPInter_interactions_file_path
        else:
            file_path = self.NPInter_interactions_old_file_path

        table = pd.read_table(file_path,
                              usecols=["ncType", "ncIdentifier", "ncName", "prType", "prIdentifier",
                                       "InteractionPartner", "organism", "tag", "interClass", "interLevel"])
        table = table[table["organism"] == "Homo sapiens"]
        table = table[table["interLevel"] == "RNA-RNA"]
        # table = table[table["interClass"].isin(["binding;regulatory", "regulatory"])]
        table["InteractionPartner"] = table["InteractionPartner"].str.replace("-3p.*|-5p.*", "")
        table["InteractionPartner"] = table["InteractionPartner"].str.replace("hsa-miR", "hsa-mir")
        table["InteractionPartner"] = table["InteractionPartner"].str.replace("miR", "hsa-mir")

        NPInter_ncRNA_RNA_regulatory_network = nx.from_pandas_edgelist(table, source='ncName',
                                                                       target='InteractionPartner',
                                                                       # edge_attr=["tag", "interClass"],
                                                                       create_using=nx.DiGraph())
        if rename_dict is not None:
            NPInter_ncRNA_RNA_regulatory_network = nx.relabel_nodes(NPInter_ncRNA_RNA_regulatory_network, rename_dict)
        return NPInter_ncRNA_RNA_regulatory_network.edges(data=data)

    def process_lncRNome_gene_info(self, lncRNome_folder_path):
        self.lnRNome_genes_info_path = os.path.join(lncRNome_folder_path, "general_information.txt")

        self.lnRNome_genes_info = pd.read_table(self.lnRNome_genes_info_path, header=0,
                                                usecols=["Gene Name", "Transcript Name", "Transcript Type", "Location", "Strand"])

    def process_lncrnadisease_associations(self, lncrnadisease_folder_path):
        self.lncrnadisease_folder_path = lncrnadisease_folder_path
        self.lncrnadisease_associations_path = os.path.join(lncrnadisease_folder_path,
                                                            "lncRNA-disease_association_v2017.txt")

        self.lncrnadisease_info = pd.read_table(self.lncrnadisease_associations_path, header=None, sep="\t")
        self.lncrnadisease_info.columns = ["LncRNA name", "Disease name", "Dysfunction type", "Description", "Chr",
                                           "Start", "End", "Strand", "Species", "Alias", "Sequence", "Reference"]
        self.lncrnadisease_info = self.lncrnadisease_info[self.lncrnadisease_info["Species"] == "Human"]

        self.lncrnadisease_info["Disease name"] = self.lncrnadisease_info["Disease name"].str.lower()

    def process_lncrna2target_interactions(self, lncrna2target_folder_path):
        self.lncrna2target_folder_path = lncrna2target_folder_path
        self.lncrna2target_high_throughput_table_path = os.path.join(self.lncrna2target_folder_path,
                                                                     "lncRNA_target_from_high_throughput_experiments.txt")
        self.lncrna2target_low_throughput_table_path = os.path.join(self.lncrna2target_folder_path,
                                                                    "lncRNA_target_from_low_throughput_experiments.xlsx")


    def get_lncrna2target_high_throughput_interactions(self, data=True, rename_dict=None):
        table = pd.read_table(self.lncrna2target_high_throughput_table_path, low_memory=True)
        table = table[table["species_id"] == 9606]
        table["lncrna_symbol"] = table["lncrna_symbol"].str.upper().replace("LINC", "")
        table["gene_symbol"] = table["gene_symbol"].str.upper()
        self.lncrna2target_high_throughput_df = table
        self.lncrna2target_high_throughput_network = nx.from_pandas_edgelist(self.lncrna2target_high_throughput_df, source='lncrna_symbol',
                                                                             target='gene_symbol',
                                                                             edge_attr=["P_Value", "direction"],
                                                                             create_using=nx.DiGraph())
        if rename_dict is not None:
            self.lncrna2target_high_throughput_network = nx.relabel_nodes(self.lncrna2target_high_throughput_network, rename_dict)
        return self.lncrna2target_high_throughput_network.edges(data=data)

    def get_lncrna2target_low_throughput_interactions(self, data=True, rename_dict=None):
        table = pd.read_excel(self.lncrna2target_low_throughput_table_path)
        table = table[table["Species"] == "Homo sapiens"]
        table["Target_official_symbol"] = table["Target_official_symbol"].str.replace("(?i)(mir)", "hsa-mir-")
        table["Target_official_symbol"] = table["Target_official_symbol"].str.replace("--", "-")
        table["Target_official_symbol"].apply(lambda x: x.lower() if "mir" in x.lower() else x.upper())
        table["GENCODE_gene_name"] = table["GENCODE_gene_name"].str.upper()
        self.lncrna2target_low_throughput_df = table
        self.lncrna2target_low_throughput_network = nx.from_pandas_edgelist(self.lncrna2target_low_throughput_df,
                                                                            source='GENCODE_gene_name',
                                                                            target='Target_official_symbol',
                                                                            create_using=nx.DiGraph())
        if rename_dict is not None:
            self.lncrna2target_low_throughput_network = nx.relabel_nodes(self.lncrna2target_low_throughput_network, rename_dict)
        return self.lncrna2target_low_throughput_network.edges(data=data)

    def process_lncRInter_interactions(self, lncRInter_folder_path):
        self.lncRInter_folder_path = lncRInter_folder_path
        self.lncRInter_table_path = os.path.join(lncRInter_folder_path, "human_interactions.txt")

    def get_lncRInter_interactions(self, data=True, rename_dict=None):
        table = pd.read_table(self.lncRInter_table_path)
        table = table[table["Organism"] == "Homo sapiens"]
        table.loc[table["Interacting partner"].str.contains("MIR"), "Interacting partner"] = table.loc[
            table["Interacting partner"].str.contains("MIR"), "Interacting partner"].str.lower()
        table["Interacting partner"] = table["Interacting partner"].str.replace("mirlet", "hsa-let-")
        table["Interacting partner"] = table["Interacting partner"].str.replace("mir", "hsa-mir-")
        table["Interacting partner"][table["Interacting partner"].str.contains(r"[mir|let]\-[\d]+[a-z]+[\d]+")] = \
            table["Interacting partner"][table["Interacting partner"].str.contains(r"[mir|let]\-[\d]+[a-z]+[\d]+")].apply(
            lambda x: x[:-1] + "-" + x[-1])
        self.lncRInter_df = table
        self.lncRInter_network = nx.from_pandas_edgelist(self.lncRInter_df, source='lncrna',
                                                             target='Interacting partner',
                                                             edge_attr=["Interaction Class", "Interaction Mode", "Tissue", "Phenotype"],
                                                             create_using=nx.DiGraph())
        if rename_dict is not None:
            self.lncRInter_network = nx.relabel_nodes(self.lncRInter_network, rename_dict)
        return self.lncRInter_network.edges(data=data)

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
        self.noncode_func_df["NONCODE Transcript ID"] = self.noncode_func_df.index.map(pd.Series(transcript2gene_df['NONCODE Transcript ID'].values,
                                                                                                 index=transcript2gene_df['NONCODE Gene ID']).to_dict())

        # Convert NONCODE transcript ID to gene names
        source_gene_names_df = source_df[source_df["name type"] == "NAME"].copy()

        self.noncode_func_df["Gene Name"] = self.noncode_func_df["NONCODE Transcript ID"].map(
            pd.Series(source_gene_names_df['Gene ID'].values,
                      index=source_gene_names_df['NONCODE Transcript ID']).to_dict())


    def process_RNAcentral_annotation_info(self, RNAcentral_folder_path):
        self.RNAcentral_annotation_file_path = os.path.join(RNAcentral_folder_path, "rnacentral_rfam_annotations.tsv")
        self.RNAcentral_gencode_id_file_path = os.path.join(RNAcentral_folder_path, "gencode.tsv")

        go_terms = pd.read_table(self.RNAcentral_annotation_file_path,
                                 low_memory=True, header=None, names=["RNAcentral id", "GO terms", "Rfams"])
        go_terms["RNAcentral id"] = go_terms["RNAcentral id"].str.split("_", expand=True)[0]

        gencode_id = pd.read_table(self.RNAcentral_gencode_id_file_path,
                                   low_memory=True, header=None,
                                   names=["RNAcentral id", "database", "external id", "species", "RNA type",
                                          "gene symbol"])
        gencode_id = gencode_id[gencode_id["species"] == 9606]

        lnc_go_terms = go_terms[go_terms["RNAcentral id"].isin(gencode_id["RNAcentral id"])].groupby("RNAcentral id")[
            "GO terms"].apply(lambda x: "|".join(x.unique()))
        lnc_rfams = go_terms[go_terms["RNAcentral id"].isin(gencode_id["RNAcentral id"])].groupby("RNAcentral id")[
            "Rfams"].apply(lambda x: "|".join(x.unique()))

        gencode_id["GO terms"] = gencode_id["RNAcentral id"].map(lnc_go_terms.to_dict())
        gencode_id["Rfams"] = gencode_id["RNAcentral id"].map(lnc_rfams.to_dict())
        gencode_id = gencode_id[gencode_id["GO terms"].notnull() | gencode_id["Rfams"].notnull()]

        self.RNAcentral_annotations = gencode_id

    def get_GENCODE_lncRNA_sequence_data(self, import_sequences, replace_U2T=True):
        seq_dict = {}
        for record in SeqIO.parse(self.GENCODE_LncRNA_sequence_file_path, "fasta"):
            # gene_id = record.id.split("|")[1]
            gene_name = record.id.split("|")[5]

            sequence_str = str(record.seq)
            if replace_U2T:
                sequence_str = sequence_str.replace("U", "T")
            if import_sequences == "shortest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) > len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif import_sequences == "longest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) < len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif import_sequences == "multi":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = [sequence_str, ]
                else:
                    seq_dict[gene_name].append(sequence_str)
            else:
                seq_dict[gene_name] = sequence_str

        return seq_dict

    def process_genes_info(self):
        # Merge lncrnadisease associations database
        self.gene_info["Disease association"] = self.gene_info["Gene Name"].map(
            self.lncrnadisease_info.groupby("LncRNA name")["Disease name"].apply('|'.join).to_dict())

        # Add RNACentral GO term and Rfam family
        self.gene_info["GO Terms"] = self.gene_info["Gene Name"].map(
            pd.Series(self.RNAcentral_annotations['GO terms'].values,
                      index=self.RNAcentral_annotations['gene symbol']).to_dict())
        self.gene_info["Rfams"] = self.gene_info["Gene Name"].map(
            pd.Series(self.RNAcentral_annotations['Rfams'].values,
                      index=self.RNAcentral_annotations['gene symbol']).to_dict())

        # Change index of genes info to gene names
        self.gene_info.index = self.get_genes_list() # Assuming the entries are ordered correctly
        self.features = list(OrderedDict.fromkeys(self.features))
        self.gene_info = self.gene_info[~self.gene_info.index.duplicated(keep='first')] # Remove duplicate genes


    def get_genes_info(self):
        return self.gene_info



class GeneExpression(ExpressionData):
    def __init__(self, cohort_name, file_path, columns="GeneSymbol|TCGA", key="GeneSymbol",
                 log2_transform=True, import_sequences="longest", replace_U2T=True):
        super().__init__(cohort_name, file_path, columns=columns, key=key, import_sequences=import_sequences, replace_U2T=replace_U2T,
                         log2_transform=log2_transform)

    def process_GENCODE_transcript_data(self, gencode_folder_path):
        self.GENCODE_transcript_sequence_file_path = os.path.join(gencode_folder_path, "gencode.v29.transcripts.fa")

    def process_targetScan_gene_info(self, targetScan_gene_info_path, human_only=True):
        self.targetScan_gene_info_path = targetScan_gene_info_path
        self.targetScan_genes_info = pd.read_table(self.targetScan_gene_info_path,
                                                   usecols=["Transcript ID", "Gene ID", "Species ID", "Gene symbol",
                                                            "Gene description", "3P-seq tags"])

        self.targetScan_genes_info["Gene description"] = self.targetScan_genes_info["Gene description"].str.replace(" \[.*\]","")

        if human_only:
            self.targetScan_genes_info = self.targetScan_genes_info[self.targetScan_genes_info["Species ID"] == 9606]
        self.targetScan_genes_info.drop(columns=["Species ID"], inplace=True)


    def process_HUGO_protein_coding_genes_info(self, hugo_protein_gene_names_path):
        self.hugo_protein_gene_names_path = hugo_protein_gene_names_path
        self.hugo_protein_genes_info = pd.read_table(self.hugo_protein_gene_names_path,
                                                     usecols=["symbol", "locus_type", "gene_family", "gene_family_id",
                                                              "location"])


    def process_GO_genes_info(self, gene_ontology_folder_path):
        self.gene_ontology_file_path = os.path.join(gene_ontology_folder_path, "goa_human.gaf")


    def get_GO_genes_info(self):
        lines = []
        with open(self.gene_ontology_file_path) as file:
            l = GOA.gafiterator(file)
            for line in l:
                lines.append(line)
        go_df = pd.DataFrame(lines)
        return go_df

    def process_RegNet_gene_regulatory_network(self, grn_file_path):
        self.regnet_grn_file_path = grn_file_path

    def process_biogrid_GRN_edgelist(self, biogrid_folder_path):
        self.biogrid_interactions_path = os.path.join(biogrid_folder_path, "BIOGRID-ALL-3.4.162.tab2.txt")

    def get_GENCODE_transcript_data(self, import_sequences, replace_U2T):
        seq_dict = {}
        locus_type_dict = {}
        for record in SeqIO.parse(self.GENCODE_transcript_sequence_file_path, "fasta"):
            gene_name = record.id.split("|")[5]
            sequence_str = str(record.seq)
            if replace_U2T:
                sequence_str = sequence_str.replace("U", "T")

            if import_sequences == "shortest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) > len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif import_sequences == "longest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) < len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif import_sequences == "multi":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = [sequence_str, ]
                else:
                    seq_dict[gene_name].append(sequence_str)
            else:
                seq_dict[gene_name] = sequence_str

            # add locus type
            if ~(gene_name in locus_type_dict):
                locus_type_dict[gene_name] = record.id.split("|")[7]
            else:
                locus_type_dict[gene_name] = locus_type_dict[gene_name] + "|" + record.id.split("|")[7]

        return seq_dict, locus_type_dict

    def process_DisGeNET_gene_disease_associations(self, disgenet_folder_path):
        self.disgenet_folder_path = disgenet_folder_path
        self.disgenet_curated_gene_disease_file_path = os.path.join(disgenet_folder_path,
                                                                    "curated_gene_disease_associations.tsv")
        self.disgenet_all_gene_disease_file_path = os.path.join(disgenet_folder_path,
                                                                "all_gene_disease_associations.tsv")

        self.disgenet_curated_gene_disease = pd.read_table(self.disgenet_curated_gene_disease_file_path,
                                                           usecols=["geneSymbol", "diseaseName", "score"])
        self.disgenet_all_gene_disease = pd.read_table(self.disgenet_all_gene_disease_file_path,
                                                       usecols=["geneSymbol", "diseaseName", "score"])

        self.disgenet_curated_gene_disease["diseaseName"] = self.disgenet_curated_gene_disease[
            "diseaseName"].str.lower()
        self.disgenet_all_gene_disease["diseaseName"] = self.disgenet_all_gene_disease["diseaseName"].str.lower()

    def process_starBase_RNA_RNA_interactions(self, starbase_folder_path):
        self.starbase_rna_rna_interaction_table_path = os.path.join(starbase_folder_path, "starbase_3.0_rna_rna_interactions.csv")

    def process_genemania_interactions(self, genemania_folder_path):
        self.genemania_interaction_table_path = os.path.join(genemania_folder_path,
                                                              "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt")

        self.genemania_identifier_mapping_path = os.path.join(genemania_folder_path,
                                                              "identifier_mappings.txt")

    def get_genemania_RNA_RNA_interactions(self, data=True):
        interactions = pd.read_table(self.genemania_interaction_table_path, low_memory=True)
        identifier = pd.read_table(self.genemania_identifier_mapping_path)

        # Rename ENSG ID's to gene names
        identifier = identifier[identifier["Source"] == "Gene Name"]
        identifier_map = pd.Series(identifier["Name"].values, index=identifier["Preferred_Name"]).to_dict()
        interactions.replace(identifier_map, inplace=True)

        self.genemania_RNA_RNA_network = nx.from_pandas_edgelist(interactions, source='Gene_A', target='Gene_B',
                                                                edge_attr=["Weight"],
                                                                create_using=nx.DiGraph())
        return self.genemania_RNA_RNA_network.edges(data=data)



    def get_starBase_RNA_RNA_interactions(self, data=True, min_interactionNum=1, min_expNum=1):
        df = pd.read_csv(self.starbase_rna_rna_interaction_table_path, header=0)

        df.loc[df["pairGeneType"]=="miRNA", "pairGeneName"] = df[df["pairGeneType"]=="miRNA"][
            "pairGeneName"].str.lower()
        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[df["pairGeneType"] == "miRNA"][
            "pairGeneName"].str.replace("-3p.*|-5p.*", "")
        df = df[df["interactionNum"] >= min_interactionNum]
        df = df[df["expNum"] >= min_expNum]

        self.starBase_RNA_RNA_network = nx.from_pandas_edgelist(df, source='geneName', target='pairGeneName',
                                                                edge_attr=["interactionNum"],
                                                                     create_using=nx.DiGraph())
        return self.starBase_RNA_RNA_network.edges(data=data)

    def get_RegNet_GRN_edgelist(self, regnet_grn_file_path=None):
        if regnet_grn_file_path is not None:
            self.regnet_grn_file_path = regnet_grn_file_path

        grn_df = pd.read_table(self.regnet_grn_file_path, header=None)

        # Since RegNet GRN contains miRNA and TF regulatory interactions
        # hsa-miR-* microRNA gene names will be mapped to hsa-mir-*
        grn_df[0] = grn_df[0].map(lambda x: x.lower() if ("hsa-miR" in x) else x)

        regnet_grn_network = nx.from_pandas_edgelist(grn_df, source=0, target=2, create_using=nx.DiGraph())

        return regnet_grn_network.edges(data=True)

    def get_BioGRID_GRN_edgelist(self, data=True, biogrid_interactions_file_path=None, rename_dict=None):
        if biogrid_interactions_file_path is not None:
            self.biogrid_interactions_path = biogrid_interactions_file_path

        biogrid_df = pd.read_table(self.biogrid_interactions_path, na_values=["-"],
                                   usecols=['Official Symbol Interactor A',
                                            'Official Symbol Interactor B', 'Organism Interactor A', 'Score',
                                            'Throughput', 'Qualifications', 'Modification', 'Phenotypes'],
                                   low_memory=True)

        biogrid_df = biogrid_df[biogrid_df["Organism Interactor A"] == 9606]
        # biogrid_df = biogrid_df[biogrid_df["Throughput"] == "High Throughput"]

        biogrid_grn = nx.from_pandas_edgelist(biogrid_df, source='Official Symbol Interactor A',
                                                   target='Official Symbol Interactor B', create_using=nx.DiGraph())
        if rename_dict is not None:
            biogrid_grn = nx.relabel_nodes(biogrid_grn, rename_dict)
        return biogrid_grn.edges(data=data) # TODO add biogrid GRN edge data?

    def process_genes_info(self, curated_gene_disease_assocs_only=True):
        self.gene_info = pd.DataFrame(index=self.get_genes_list())

        self.gene_info.index.name = "Gene symbol"
        self.gene_info = self.gene_info.join(self.targetScan_genes_info.groupby("Gene symbol").first(), on="Gene symbol", how="left")

        self.gene_info.index.name = "symbol"
        self.gene_info = self.gene_info.join(self.hugo_protein_genes_info.groupby("symbol").first(), on="symbol", how="left")

        transcript_seq, gene_type = self.get_GENCODE_transcript_data(self.import_sequences, self.replace_U2T)
        self.gene_info["locus_type"] = self.gene_info.index.map(gene_type)
        self.gene_info["Transcript sequence"] = self.gene_info.index.map(transcript_seq)

        go_df = self.get_GO_genes_info()
        self.gene_info["GO Terms"] = self.gene_info.index.map(go_df.groupby("DB_Object_Symbol")["GO_ID"].apply(
            lambda x: "|".join(x.unique())).to_dict())

        if curated_gene_disease_assocs_only:
            self.gene_info["Disease association"] = self.gene_info.index.map(
                self.disgenet_curated_gene_disease.groupby("geneSymbol")["diseaseName"].apply('|'.join).to_dict())
        else:
            self.gene_info["Disease association"] = self.gene_info.index.map(
                self.disgenet_all_gene_disease.groupby("geneSymbol")["diseaseName"].apply('|'.join).to_dict())

        # Process gene location info
        self.gene_info["Chromosome"] = "chr" + self.gene_info["location"].str.split("p|q", expand=True)[0]
        self.gene_info["Chromosome arm"] = self.gene_info["location"].str.extract(r'(?P<arm>[pq])', expand=True)
        self.gene_info["Chromosome region"] = self.gene_info["location"].str.split("[pq.-]", expand=True)[1]
        self.gene_info["Chromosome band"] = self.gene_info["location"].str.split("[pq.-]", expand=True)[2]

        self.gene_info["Transcript length"] = self.gene_info["Transcript sequence"].apply(
            lambda x: len(x) if type(x) is str else None)

        self.gene_info["3P-seq tags"] = self.gene_info["3P-seq tags"].astype("O")

    def get_genes_info(self):
        return self.gene_info


class MiRNAExpression(ExpressionData):
    def __init__(self, cohort_name, file_path, columns="GeneSymbol|TCGA", key="GeneSymbol", log2_transform=True, import_sequences="longest", replace_U2T=True):
        super().__init__(cohort_name, file_path, columns=columns, key=key, import_sequences=import_sequences, replace_U2T=replace_U2T,
                         log2_transform=log2_transform)

    def process_mirbase_data(self, mirbase_folder_path):
        self.mirbase_aliases_file_path = os.path.join(mirbase_folder_path, "aliases.txt")
        self.mirbase_mir_seq_file_path = os.path.join(mirbase_folder_path, "hairpin.fa")

    def get_mirbase_hairpin_sequence_data(self, import_sequences, replace_U2T=True):
        seq_dict = {}
        for record in SeqIO.parse(self.mirbase_mir_seq_file_path, "fasta"):
            gene_name = str(record.id)
            sequence_str = str(record.seq)
            if replace_U2T:
                sequence_str = sequence_str.replace("U", "T")

            if import_sequences == "shortest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) > len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif import_sequences == "longest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) < len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif import_sequences == "multi":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = [sequence_str, ]
                else:
                    seq_dict[gene_name].append(sequence_str)
            else:
                seq_dict[gene_name] = sequence_str

        return seq_dict

    def process_target_scan(self, targetScan_folder_path):
        self.targetScan_miR_family_info_path = os.path.join(targetScan_folder_path,"miR_Family_Info.txt")
        self.targetScan_predicted_targets_path = os.path.join(targetScan_folder_path, "Predicted_Targets_Info.default_predictions.txt")
        self.targetScan_predicted_targets_context_score_path = os.path.join(targetScan_folder_path, "Predicted_Targets_Context_Scores.default_predictions.txt")

        self.process_targetscan_mirna_family()
        self.process_targetScan_mirna_target_interactions()
        self.process_targetScan_mirna_target_interactions_context_score()

    def process_targetscan_mirna_family(self, human_only=True):
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
        self.targetScan_family_df.drop('MiRBase ID', axis=1, inplace=True)

    def process_targetScan_mirna_target_interactions(self):
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
        self.targetScan_df = targetScan_df

    def process_targetScan_mirna_target_interactions_context_score(self):
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
        self.targetScan_context_df = targetScan_context_df

    def process_mirnadisease_associations(self, HMDD_miRNAdisease_path):
        self.HMDD_miRNAdisease_path = HMDD_miRNAdisease_path
        self.mirnadisease_association_path = os.path.join(HMDD_miRNAdisease_path, "miRNA_disease.txt")

        self.mirnadisease = pd.read_table(self.mirnadisease_association_path, header=None, sep="\t")
        self.mirnadisease.columns = ["index", "miRNA name", "Disease name", "Reference", "Description"]

        self.mirnadisease["Disease name"] = self.mirnadisease["Disease name"].str.lower()

    def process_miRTarBase_miRNA_target_interactions(self, miRTarBase_path):
        self.miRTarBase_path = miRTarBase_path
        self.miRTarBase_MTI_path = os.path.join(self.miRTarBase_path, "miRTarBase_MTI.xlsx")
        self.miRTarBase_MTI_old_path = os.path.join(self.miRTarBase_path, "miRTarBase_MTI_v6.xlsx")

    def process_HUGO_miRNA_gene_info(self, HUGO_folder_path):
        self.HUGO_miRNA_gene_info_path = os.path.join(HUGO_folder_path, "RNA_micro.txt")

        df = pd.read_table(self.HUGO_miRNA_gene_info_path,
                           usecols=["alias_symbol", "ensembl_gene_id", "location"])
        df = df[df["alias_symbol"].notna()]

        df_1 = pd.concat([Series(data=row["location"], index=row['alias_symbol'].split('|'))
                          for _, row in df.iterrows()]).reset_index()
        df_2 = pd.concat([Series(data=row["ensembl_gene_id"], index=row['alias_symbol'].split('|'))
                          for _, row in df.iterrows()]).reset_index()

        df_1.columns = ["MiRBase ID", "location"]
        df_2.columns = ["MiRBase ID", "ensembl_gene_id"]

        self.HUGO_miRNA_gene_info_df = pd.merge(df_1, df_2, how="inner", on="MiRBase ID")
        self.HUGO_miRNA_gene_info_df.index = self.HUGO_miRNA_gene_info_df["MiRBase ID"]

    def process_RNAcentral_annotation_info(self, RNAcentral_folder_path):
        self.RNAcentral_annotation_file_path = os.path.join(RNAcentral_folder_path, "rnacentral_rfam_annotations.tsv")
        self.RNAcentral_mirbase_id_file_path = os.path.join(RNAcentral_folder_path, "mirbase.tsv")

        # Add tables for ID look up
        mirbase_id = pd.read_table(self.RNAcentral_mirbase_id_file_path,
                                   low_memory=True, header=None,
                                   names=["RNAcentral id", "database", "mirbase id", "species", "RNA type",
                                          "gene name"],
                                   index_col="mirbase id")
        mirbase_id = mirbase_id[mirbase_id["species"] == 9606]

        mirbase_name = pd.read_table(self.mirbase_aliases_file_path,
                                     low_memory=True, header=None, names=["mirbase id", "miRNA name"], dtype="O")
        mirbase_name = mirbase_name.join(mirbase_id, on="mirbase id", how="inner")


        # filter GO terms and Rfam annotation for miRNAs
        go_terms = pd.read_table(self.RNAcentral_annotation_file_path,
                                 low_memory=True, header=None, names=["RNAcentral id", "GO terms", "Rfams"])
        go_terms["RNAcentral id"] = go_terms["RNAcentral id"].str.split("_", expand=True)[0]

        mir_go_terms = go_terms[go_terms["RNAcentral id"].isin(mirbase_name["RNAcentral id"])].groupby("RNAcentral id")[
            "GO terms"].apply(lambda x: "|".join(x.unique()))
        mir_rfam = go_terms[go_terms["RNAcentral id"].isin(mirbase_name["RNAcentral id"])].groupby("RNAcentral id")[
            "Rfams"].apply(lambda x: "|".join(x.unique()))

        mirbase_name["GO terms"] = mirbase_name["RNAcentral id"].map(mir_go_terms.to_dict())
        mirbase_name["Rfams"] = mirbase_name["RNAcentral id"].map(mir_rfam.to_dict())

        # Expanding miRNA names in each MirBase Ascension ID
        s = mirbase_name.apply(lambda x: pd.Series(x['miRNA name'].split(";")[:-1]), axis=1).stack().reset_index(
            level=1, drop=True)
        s.name = "miRNA name"
        mirbase_name = mirbase_name.drop('miRNA name', axis=1).join(s)
        mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        self.RNAcentral_annotations = mirbase_name


    def get_targetScan_miRNA_target_interaction(self):
        if self.targetScan_df is None:
            raise Exception("must first run process_target_scan(mirna_list, gene_symbols)")

        mir_target_network = nx.from_pandas_edgelist(self.targetScan_df, source="MiRBase ID", target="Gene Symbol",
                                                     create_using=nx.DiGraph())
        return mir_target_network.edges(data=True)

    def get_targetScan_miRNA_target_interactions_context_edgelist(self, data=True):
        mirna_target_interactions = self.targetScan_context_df.copy()
        mirna_target_interactions["weighted context++ score percentile"] = \
            mirna_target_interactions["weighted context++ score percentile"].apply(func=lambda x: x / 100.0)
        mirna_target_interactions.rename(columns={"weighted context++ score percentile": "weight",
                                                  "MiRBase ID": "MIR",
                                                  "Gene Symbol": "GE"}, inplace=True)
        mir_target_network = nx.from_pandas_edgelist(mirna_target_interactions, source="MIR", target="GE",
                                                     edge_attr="weight", create_using=nx.DiGraph())
        return mir_target_network.edges(data=data)

    def get_miRTarBase_miRNA_target_interaction(self, use_latest=True, data=True, rename_dict=None):
        if use_latest:
            table = pd.read_excel(self.miRTarBase_MTI_path)
        else:
            table = pd.read_excel(self.miRTarBase_MTI_old_path)

        table = table[table["Species (Target Gene)"] == "Homo sapiens"]

        table['miRNA'] = table['miRNA'].str.lower()
        table['miRNA'] = table['miRNA'].str.replace("-3p.*|-5p.*", "")
        self.miRTarBase_df = table

        if self.miRTarBase_df is None:
            raise Exception("must first run process_miRTarBase_miRNA_target_interactions")

        mir_target_network = nx.from_pandas_edgelist(self.miRTarBase_df, source="miRNA", target="Target Gene",
                                                     edge_attr=["Support Type"],
                                                     create_using=nx.DiGraph())
        if rename_dict is not None:
            mir_target_network = nx.relabel_nodes(mir_target_network, rename_dict)
        return mir_target_network.edges(data=data)

    def process_genes_info(self):
        self.gene_info = pd.DataFrame(index=self.get_genes_list())
        self.gene_info.index.name = "MiRBase ID"

        self.gene_info = self.gene_info.join(self.targetScan_family_df.groupby("MiRBase ID").first(), on="MiRBase ID", how="left")
        self.gene_info = self.gene_info.join(self.HUGO_miRNA_gene_info_df, on="MiRBase ID", how="left")

        self.gene_info["Disease association"] = self.gene_info.index.map(
            self.mirnadisease.groupby("miRNA name")["Disease name"].apply('|'.join).to_dict())

        self.gene_info["locus_type"] = "microRNA"

        # self.gene_info.rename(columns={'Mature sequence': 'Transcript sequence'}, inplace=True)
        self.gene_info["Transcript sequence"] = self.gene_info["MiRBase ID"].map(
            self.get_mirbase_hairpin_sequence_data(self.import_sequences, self.replace_U2T))
        self.gene_info["Transcript length"] = self.gene_info["Transcript sequence"].apply(
            lambda x: len(x) if type(x) is str else None)

        # Process gene location info
        self.gene_info["Chromosome"] = "chr" + self.gene_info["location"].str.split("p|q", expand=True)[0]
        self.gene_info["Chromosome arm"] = self.gene_info["location"].str.extract(r'(?P<arm>[pq])', expand=True)
        self.gene_info["Chromosome region"] = self.gene_info["location"].str.split("[pq.-]", expand=True)[1]
        self.gene_info["Chromosome band"] = self.gene_info["location"].str.split("[pq.-]", expand=True)[2]

        self.gene_info["Transcript length"] = self.gene_info["Transcript sequence"].apply(
            lambda x: len(x) if type(x) is str else None)

        # Process Annotation data
        self.gene_info["GO Terms"] = self.gene_info["MiRBase ID"].map(
            pd.Series(self.RNAcentral_annotations['GO terms'].values,
                      index=self.RNAcentral_annotations['miRNA name']).to_dict())
        self.gene_info["Rfams"] = self.gene_info["MiRBase ID"].map(
            pd.Series(self.RNAcentral_annotations['Rfams'].values,
                      index=self.RNAcentral_annotations['miRNA name']).to_dict())

    def get_genes_info(self):
        return self.gene_info


class ProteinExpression(ExpressionData):
    def __init__(self, cohort_name, file_path, columns="GeneSymbol|TCGA", key="GeneSymbol", import_sequences="longest", log2_transform=True):
        super().__init__(cohort_name, file_path, columns=columns, key=key, import_sequences=import_sequences, log2_transform=log2_transform)

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_edgelist(HPRD_PPI, source=0, target=3,
                                          create_using=nx.DiGraph())

    def get_HPRD_PPI_network_edgelist(self):
        return self.HPRD_PPI_network.edges(data=True)

    def process_STRING_PPI_network(self, ppi_data_file_path):
        pass


