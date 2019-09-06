from collections import OrderedDict

import networkx as nx
import numpy as np
from Bio.UniProt import GOA
from pandas import Series
import pandas as pd
import dask.dataframe as dd

from openomics.database.annotation import *


class ExpressionData(object):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        """
        .. class:: ExpressionData
        An abstract class that handles importing of any quantitative -omics data that is in a table format (e.g. csv, tsv, excel). Pandas will load the DataFrame from file with the user-specified columns and genes column name, then tranpose it such that the rows are samples and columns are gene/transcript/peptides.
        The user will also specify the index argument, which specifies if the genes are ensembl genes ID or gene name, or transcripts id/names. The user should be careful about choosing the right genes index which makes it easier to annotate functional, sequence, and interaction data to it.
        The dataframe should only contain numeric values besides the genes_col_name and the sample barcode id indices.
        Args:
            cohort_name (str): the unique cohort code name string
            index (str): {"gene_id", "transcript_id", "peptide_id", "gene_name", "trascript_name", "peptide_name"}
                Chooses the level of the gene/transcript/peptide of the genes list in this expression data. The expression DataFrame's index will be renamed to this.
            file_path (str):
                Path of the table file to import.
            columns (str): a regex string
                A regex string to import column names from the table. Columns names imported are string match, separated by "|".
            genes_col_name (str):
                Index column name of expression data which is used to index the genes list
            transposed (bool): default True
                If True, perform preprocessing steps for the table data obtained from TCGA-Assembler tool. If False, import a pandas table as-is with bcr_sample_barcode for row index, and gene names as columns
            log2_transform (bool): default False
                Whether to log2 transform the expression values
        """
        self.cohort_name = cohort_name

        if os.path.isfile(file_path) and os.path.exists(file_path):
            table = pd.read_table(file_path)
        else:
            raise FileNotFoundError(file_path)

        self.index = index
        self.expressions = self.preprocess_table(table, columns, genes_col_name, transposed)
        self.expressions.index.name = self.index

        if log2_transform:
            self.expressions = self.expressions.applymap(self.log2_transform)

        # Save samples and features for this omics data
        self.samples = self.expressions.index
        self.features = self.expressions.columns.tolist()

    def preprocess_table(self, df, columns, key, transposed):
        # type: (pd.DataFrame, str, str, bool) -> None
        """
        This function preprocesses the expression table files where columns are samples and rows are gene/transcripts
        :param transposed:
        :return:
        """
        # Filter columns
        df = df.filter(regex=columns)

        # Cut TCGA column names to sample barcode, discarding aliquot info
        df = df.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x, inplace=False)

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

        # Drop NA geneID rows
        df.dropna(axis=0, inplace=True)

        # Remove entries with unknown geneID
        df = df[df[key] != '?']

        # Transpose dataframe to patient rows and geneID columns
        df.index = df[key]
        df.drop([key], axis=1, inplace=True)

        if transposed:
            df = df.T

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

        return df

    def log2_transform(self, x):
        return np.log2(x + 1)

    def drop_genes(self, genes_to_drop):
        self.expressions.drop(genes_to_drop, axis=1, inplace=True)
        for gene in genes_to_drop:
            self.features.remove(gene)

    @classmethod
    def name(cls):
        raise NotImplementedError

    def get_genes_list(self):
        return self.features

    def get_samples_list(self):
        return self.samples


class LncRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        """
        :param file_path: Path to the lncRNA expression data, downloaded from http://ibl.mdanderson.org/tanric/_design/basic/index.html

        Args:
            level:
        """
        super(LncRNA, self).__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                                     transposed=transposed, log2_transform=log2_transform)

    @classmethod
    def name(cls):
        return cls.__name__

    def preprocess_table(self, df, columns, key, transposed):
        # type: (pd.DataFrame, str, str, bool) -> None
        """
        Preprocess LNCRNA expression file obtained from TANRIC MDAnderson, and replace ENSEMBL gene ID to HUGO gene names (HGNC). This function overwrites the GenomicData.process_expression_table() function which processes TCGA-Assembler data.

        TANRIC LNCRNA expression values are log2 transformed
        :param transposed:
        """

        # Replacing ENSG Gene ID to the lncRNA gene symbol name
        df[key] = df[key].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
        df = df[~df[key].duplicated(keep='first')] # Remove duplicate genes

        # Drop NA gene rows
        df.dropna(axis=0, inplace=True)

        # Transpose matrix to patients rows and genes columns
        df.index = df[key]
        df = df.T.iloc[1:, :]

        # Change index string to bcr_sample_barcode standard
        def change_patient_barcode(s):
            if "Normal" in s:
                return s[s.find('TCGA'):] + "-11A"
            elif "Tumor" in s:
                return s[s.find('TCGA'):] + "-01A"
            else:
                return s

        df.index = df.index.map(change_patient_barcode)
        df.index.name = "gene_id"

        return df

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

    # def get_lncBase_miRNA_lncRNA_predicted_interactions_edgelist(self, data=True, rename_dict=None):
    #     records = []
    #     for record in SeqIO.parse(
    #             self.lncBase_predicted_interactions_file_path,
    #             "fasta"):
    #         records.append(record.id.split(","))
    #     lncbase_df = pd.DataFrame(records, columns=["Transcript_ID", "geneId", "mirna", "miTG-score"])
    #
    #     lncbase_df["miTG-score"] = lncbase_df["miTG-score"].astype(float)
    #     lncbase_df = lncbase_df[lncbase_df["miTG-score"] >= 0.9]
    #
    #     lncbase_df["geneId"] = lncbase_df["geneId"].str.replace("(\(.+|\))", "")
    #     lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("(\(.+|\))", "")
    #     lncbase_df["mirna"] = lncbase_df["mirna"].str.lower()
    #     lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("-3p.*|-5p.*", "")
    #
    #     lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source='mirna', target='geneId',
    #                                                            edge_attr=["miTG-score"],
    #                                                            create_using=nx.DiGraph())
    #     if rename_dict is not None:
    #         lncBase_lncRNA_miRNA_network = nx.relabel_nodes(lncBase_lncRNA_miRNA_network, rename_dict)
    #
    #     return lncBase_lncRNA_miRNA_network.edges(data=data)

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



    def process_genes_info(self):
        # Merge lncrnadisease associations database
        self.annotations["Disease association"] = self.annotations["Gene Name"].map(
            self.lncrnadisease_info.groupby("LncRNA name")["Disease name"].apply('|'.join).to_dict())

        # Change index of genes info to gene names
        self.annotations = self.annotations[~self.annotations.index.duplicated(keep='first')] # Remove duplicate genes


    def get_annotations(self):
        return self.annotations



class MessengerRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super(MessengerRNA, self).__init__(cohort_name, index, file_path, columns=columns,
                                           genes_col_name=genes_col_name,
                                           transposed=transposed, log2_transform=log2_transform)

    @classmethod
    def name(cls):
        return cls.__name__

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


    def process_genes_info(self, curated_gene_disease_assocs_only=True):
        self.gene_info = pd.DataFrame(index=self.get_genes_list())

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

    def get_annotations(self):
        return self.gene_info


class MicroRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super(MicroRNA, self).__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                                       transposed=transposed, log2_transform=log2_transform)

    @classmethod
    def name(cls):
        return cls.__name__

    def process_mirnadisease_associations(self, HMDD_miRNAdisease_path):
        self.HMDD_miRNAdisease_path = HMDD_miRNAdisease_path
        self.mirnadisease_association_path = os.path.join(HMDD_miRNAdisease_path, "miRNA_disease.txt")

        self.mirnadisease = pd.read_table(self.mirnadisease_association_path, header=None, sep="\t")
        self.mirnadisease.columns = ["index", "miRNA name", "Disease name", "Reference", "Description"]

        self.mirnadisease["Disease name"] = self.mirnadisease["Disease name"].str.lower()

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


    def process_genes_info(self):
        self.gene_info = pd.DataFrame(index=self.get_genes_list())
        self.gene_info.index.name = "MiRBase ID"

        self.gene_info = self.gene_info.join(self.targetScan_family_df.groupby("MiRBase ID").first(), on="MiRBase ID", how="left")
        self.gene_info = self.gene_info.join(self.HUGO_miRNA_gene_info_df, on="MiRBase ID", how="left")

        self.gene_info["Disease association"] = self.gene_info.index.map(
            self.mirnadisease.groupby("miRNA name")["Disease name"].apply('|'.join).to_dict())


    def get_annotations(self):
        return self.gene_info


