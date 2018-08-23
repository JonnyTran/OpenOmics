import os
from collections import OrderedDict

import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO
from pandas import Series

from TCGAMultiOmics.utils import GTF


class GenomicData:
    def __init__(self, cancer_type, file_path, columns="GeneSymbol|TCGA",
                 import_from_TCGA_Assembler=True, log2_transform=True):
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

    def drop_genes(self, genes_to_drop):
        self.data.drop(genes_to_drop, axis=1, inplace=True)
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
            return self.network.edges()
        else:
            print(self.__class__.__str__(), "does not have network interaction data yet. (at self.network)")
            return None


class LncRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, HGNC_lncRNA_names_file_path, GENCODE_folder_path):
        """
        :param folder_path: Path to the lncRNA expression data, downloaded from http://ibl.mdanderson.org/tanric/_design/basic/index.html
        :param lncrna_names_file_path: Path to the HGNC_RNA_long_non-coding.txt file to map ensembl gene id to lncRNA names
        """
        file_path = os.path.join(folder_path, "TCGA-rnaexpr.tsv")
        self.HGNC_lncRNA_names_path = HGNC_lncRNA_names_file_path
        self.GENCODE_LncRNA_gtf_file_path = os.path.join(GENCODE_folder_path, "gencode.v28.long_noncoding_RNAs.gtf")
        self.GENCODE_LncRNA_sequence_file_path = os.path.join(GENCODE_folder_path, "gencode.v28.lncRNA_transcripts.fa")
        super().__init__(cancer_type, file_path, log2_transform=False)


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
            HGNC_lncrna_info = pd.read_table(self.HGNC_lncRNA_names_path, delimiter="\t", usecols=['symbol', 'locus_type', 'ensembl_gene_id', 'name', 'location'])
            self.HGNC_lncrna_info = HGNC_lncrna_info
            self.HGNC_lncrna_info.index = self.HGNC_lncrna_info["symbol"]
        except Exception:
            raise FileNotFoundError("Needs the file RNA_long_non-coding.txt at directory external_data/HUGO_Gene_names to process lncRNA gene info")


        # Replacing ENSG Gene ID to the lncRNA gene symbol name
        ensembl_id_to_gene_name, ensembl_id_to_transcript_id = self.get_GENCODE_lncRNA_gene_name_dict()
        hgnc_lncrna_dict = self.get_HUGO_lncRNA_gene_name_dict()
        lncrna_exp['Gene_ID'] = lncrna_exp['Gene_ID'].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
        lncrna_exp = lncrna_exp[~lncrna_exp['Gene_ID'].duplicated(keep='first')] # Remove duplicate genes


        # Preprocess genes info
        self.preprocess_genes_info(lncrna_exp['Gene_ID'], ensembl_id_to_gene_name, ensembl_id_to_transcript_id, hgnc_lncrna_dict)



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

        self.starBase_miRNA_lncRNA_network = nx.from_pandas_edgelist(grn_df, source='geneName', target='name', create_using=nx.DiGraph())

    def get_starBase_lncRNA_miRNA_interactions_edgelist(self):
        return self.starBase_miRNA_lncRNA_network.edges()

    def process_lncRNome_miRNA_binding_sites(self, lncRNome_folder_path):
        self.lnRNome_miRNA_binding_sites_path = os.path.join(lncRNome_folder_path, "miRNA_binding_sites.txt")

        df = pd.read_table(self.lnRNome_miRNA_binding_sites_path, header=0)

        df['Binding miRNAs'] = df['Binding miRNAs'].str.lower()
        df['Binding miRNAs'] = df['Binding miRNAs'].str.replace("-3p.*|-5p.*", "")

        self.lncRNome_miRNA_binding_sites_network = nx.from_pandas_edgelist(df, source='Gene Name', target='Binding miRNAs', create_using=nx.DiGraph())

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

    def process_lncrnadisease_associations(self, lncrnadisease_folder_path):
        self.lncrnadisease_folder_path = lncrnadisease_folder_path
        self.lncrnadisease_associations_path = os.path.join(lncrnadisease_folder_path,
                                                            "lncRNA-disease_association_v2017.txt")

        self.lncrnadisease_info = pd.read_table(self.lncrnadisease_associations_path, header=None, sep="\t")
        self.lncrnadisease_info.columns = ["LncRNA name", "Disease name", "Dysfunction type", "Description", "Chr",
                                           "Start", "End", "Strand", "Species", "Alias", "Sequence", "Reference"]
        self.lncrnadisease_info = self.lncrnadisease_info[self.lncrnadisease_info["Species"] == "Human"]

        self.lncrnadisease_info["Disease name"] = self.lncrnadisease_info["Disease name"].str.lower()



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
        source_gene_names_df = source_df[source_df["name type"] == "NAME"].copy()

        self.noncode_func_df["Gene Name"] = self.noncode_func_df["NONCODE Transcript ID"].map(
            pd.Series(source_gene_names_df['Gene ID'].values,
                      index=source_gene_names_df['NONCODE Transcript ID']).to_dict())


        # TODO Convert using lncpedia ID

        # # Convert NONCODE transcript ID to gencode transcript ID
        # source_gene_names_df = source_df[source_df["name type"] == "gencode"].copy()
        # self.noncode_func_df["gencode transcript id"] = self.noncode_func_df["NONCODE Transcript ID"].map(
        #     pd.Series(source_gene_names_df['Gene ID'].values,
        #               index=source_gene_names_df['NONCODE Transcript ID']).to_dict())
        #
        # # Convert NONCODE transcript ID to ensembl transcript ID
        # source_gene_names_df = source_df[source_df["name type"] == "ensembl"].copy()
        # self.noncode_func_df["ensembl transcript id"] = self.noncode_func_df["NONCODE Transcript ID"].map(
        #     pd.Series(source_gene_names_df['Gene ID'].values,
        #               index=source_gene_names_df['NONCODE Transcript ID']).to_dict())

    def get_GENCODE_lncRNA_sequence_data(self):
        lnc_seq = {}
        for record in SeqIO.parse(self.GENCODE_LncRNA_sequence_file_path, "fasta"):
            lnc_seq[record.id.split("|")[5]] = str(record.seq)

        return lnc_seq

    def process_genes_info(self):
        self.gene_info.index.name = "symbol"
        self.gene_info = self.gene_info.join(self.HGNC_lncrna_info.groupby("symbol").first(), on="symbol",
                                             how="left")

        self.gene_info.index.name = "Gene Name"
        self.gene_info = self.gene_info.join(self.lnRNome_genes_info.groupby("Gene Name").first(), on="Gene Name",
                                             how="left")

        self.gene_info.index = self.gene_info["Gene Name"]
        self.gene_info = self.gene_info.join(self.noncode_func_df.groupby("Gene Name").first(), on="Gene Name",
                                             how="left")

        # Merge transcript sequence data
        self.gene_info["Transcript sequence"] = self.gene_info["Gene Name"].map(
            self.get_GENCODE_lncRNA_sequence_data())

        # Merge lncrnadisease associations database
        self.gene_info["Disease association"] = self.gene_info["Gene Name"].map(
            self.lncrnadisease_info.groupby("LncRNA name")["Disease name"].apply('|'.join).to_dict())


        self.gene_info.index = self.get_genes_list() # Assuming the entries are ordered correctly

        self.features = list(OrderedDict.fromkeys(self.features))
        self.gene_info = self.gene_info[~self.gene_info.index.duplicated(keep='first')] # Remove duplicate genes

        # Process gene location info
        self.gene_info["Chromosome"] = self.gene_info["Location"].str.split(":", expand=True)[0]
        self.gene_info["start"] = self.gene_info["Location"].str.split(":", expand=True)[1].str.split("-", expand=True)[0]
        self.gene_info["end"] = self.gene_info["Location"].str.split(":", expand=True)[1].str.split("-", expand=True)[1]
        self.gene_info["bp length"] = self.gene_info["Location"].str.split(":", expand=True)[1].apply(
            lambda x: int(x.split("-")[1]) - int(x.split("-")[0]) if (type(x) is str) else None)

        self.gene_info["Transcript length"] = self.gene_info["Transcript sequence"].apply(lambda x: len(x) if type(x) is str else None)


        self.gene_info["start"].astype(np.float64)
        self.gene_info["end"].astype(np.float64)

        self.gene_info["locus_type"] = "RNA, long non-coding"

    def get_genes_info(self):
        return self.gene_info



class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, log2_transform=True):
        file_path = os.path.join(folder_path, "geneExp.txt")
        super().__init__(cancer_type, file_path, log2_transform=log2_transform)

    def process_GENCODE_transcript_data(self, gencode_folder_path):
        self.GENCODE_transcript_sequence_file_path = os.path.join(gencode_folder_path, "gencode.v28.transcripts.fa")

    def process_targetScan_gene_info(self, targetScan_gene_info_path, human_only=True):
        self.targetScan_gene_info_path = targetScan_gene_info_path
        self.targetScan_genes_info = pd.read_table(self.targetScan_gene_info_path, usecols=["Transcript ID", "Gene ID", "Species ID", "Gene symbol", "Gene description", "3P-seq tags"])

        self.targetScan_genes_info["Gene description"] = self.targetScan_genes_info["Gene description"].str.replace(" \[.*\]","")

        if human_only:
            self.targetScan_genes_info = self.targetScan_genes_info[self.targetScan_genes_info["Species ID"] == 9606]
        self.targetScan_genes_info.drop(columns=["Species ID"], inplace=True)


    def process_HUGO_protein_coding_genes_info(self, hugo_protein_gene_names_path):
        self.hugo_protein_gene_names_path = hugo_protein_gene_names_path
        self.hugo_protein_genes_info = pd.read_table(self.hugo_protein_gene_names_path, usecols=["symbol", "locus_type", "gene_family", "gene_family_id", "location"])


    def process_RegNet_gene_regulatory_network(self, grn_file_path):
        self.regnet_grn_file_path = grn_file_path

    def process_biogrid_GRN_edgelist(self, biogrid_folder_path):
        self.biogrid_interactions_path = os.path.join(biogrid_folder_path, "BIOGRID-ALL-3.4.162.tab2.txt")

    def get_GENCODE_transcript_data(self):
        transcript_seq = {}
        locus_type_dict = {}
        for record in SeqIO.parse(self.GENCODE_transcript_sequence_file_path, "fasta"):
            gene_name = record.id.split("|")[5]
            transcript_seq[gene_name] = str(record.seq)
            if ~(gene_name in locus_type_dict):
                locus_type_dict[gene_name] = record.id.split("|")[7]
            else:
                locus_type_dict[gene_name] = locus_type_dict[gene_name] + "|" + record.id.split("|")[7]

        return transcript_seq, locus_type_dict

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


    def get_RegNet_GRN_edgelist(self, regnet_grn_file_path=None):
        if regnet_grn_file_path is not None:
            self.regnet_grn_file_path = regnet_grn_file_path

        grn_df = pd.read_table(self.regnet_grn_file_path, header=None)

        # Since RegNet GRN contains miRNA and TF regulatory interactions
        # hsa-miR-* microRNA gene names will be mapped to hsa-mir-*
        grn_df[0] = grn_df[0].map(lambda x: x.lower() if ("hsa-miR" in x) else x)

        regnet_grn_network = nx.from_pandas_edgelist(grn_df, source=0, target=2, create_using=nx.DiGraph())

        return regnet_grn_network.edges(data=True)

    def get_BioGRID_GRN_edgelist(self, biogrid_interactions_file_path=None):
        if biogrid_interactions_file_path is not None:
            self.biogrid_interactions_path = biogrid_interactions_file_path

        biogrid_df = pd.read_table(self.biogrid_interactions_path, na_values=["-"],
                                   usecols=['Official Symbol Interactor A',
                                            'Official Symbol Interactor B', 'Organism Interactor A', 'Score',
                                            'Throughput', 'Qualifications', 'Modification', 'Phenotypes'])

        biogrid_df = biogrid_df[biogrid_df["Organism Interactor A"] == 9606]

        biogrid_grn = nx.from_pandas_edgelist(biogrid_df, source='Official Symbol Interactor A',
                                                   target='Official Symbol Interactor B', create_using=nx.DiGraph())
        return biogrid_grn.edges(data=False) # TODO add biogrid GRN edge data?

    def process_genes_info(self, curated_gene_disease_assocs_only=True):
        self.gene_info = pd.DataFrame(index=self.get_genes_list())

        self.gene_info.index.name = "Gene symbol"
        self.gene_info = self.gene_info.join(self.targetScan_genes_info.groupby("Gene symbol").first(), on="Gene symbol", how="left")

        self.gene_info.index.name = "symbol"
        self.gene_info = self.gene_info.join(self.hugo_protein_genes_info.groupby("symbol").first(), on="symbol", how="left")

        transcript_seq, gene_type = self.get_GENCODE_transcript_data()
        self.gene_info["locus_type"] = self.gene_info.index.map(gene_type)
        self.gene_info["Transcript sequence"] = self.gene_info.index.map(transcript_seq)

        if curated_gene_disease_assocs_only:
            self.gene_info["Disease association"] = self.gene_info.index.map(
                self.disgenet_curated_gene_disease.groupby("geneSymbol")["diseaseName"].apply('|'.join).to_dict())
        else:
            self.gene_info["Disease association"] = self.gene_info.index.map(
                self.disgenet_all_gene_disease.groupby("geneSymbol")["diseaseName"].apply('|'.join).to_dict())

        # Process gene location info
        self.gene_info["Chromosome"] = "Chromosome " + self.gene_info["location"].str.split("p|q", expand=True)[0]
        self.gene_info["Chromosome arm"] = self.gene_info["location"].str.extract(r'(?P<arm>[pq])', expand=True)
        self.gene_info["Chromosome region"] = self.gene_info["location"].str.split("[pq.-]", expand=True)[0]

        self.gene_info["Transcript length"] = self.gene_info["Transcript sequence"].apply(
            lambda x: len(x) if type(x) is str else None)

    def get_genes_info(self):
        return self.gene_info


class SomaticMutation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "somaticMutation_geneLevel.txt")
        super().__init__(cancer_type, file_path)


class MiRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, log2_transform=True):
        file_path = os.path.join(folder_path, "miRNAExp__RPM.txt")
        super().__init__(cancer_type, file_path, log2_transform=log2_transform)

    def process_target_scan(self, targetScan_folder_path):
        self.targetScan_miR_family_info_path = os.path.join(targetScan_folder_path,"miR_Family_Info.txt")
        self.targetScan_predicted_targets_path = os.path.join(targetScan_folder_path, "Predicted_Targets_Info.default_predictions.txt")
        self.targetScan_predicted_targets_context_score_path = os.path.join(targetScan_folder_path, "Predicted_Targets_Context_Scores.default_predictions.txt")

        self.process_targetscan_mirna_family()
        self.process_mirna_target_interactions()
        self.process_mirna_target_interactions_context_score()

    def process_targetscan_mirna_family(self, human_only=True, incremental_group_numbering=False):
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

    def process_mirna_target_interactions(self):
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

    def process_mirna_target_interactions_context_score(self):
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


    def get_miRNA_target_interaction(self):
        if self.targetScan_df is None:
            raise Exception("must first run process_target_scan(mirna_list, gene_symbols)")
        return self.targetScan_df

    def get_miRNA_target_interaction_edgelist(self):
        mirna_target_interactions = self.targetScan_context_df.copy()
        mirna_target_interactions["weighted context++ score percentile"] = \
            mirna_target_interactions["weighted context++ score percentile"].apply(func=lambda x: x / 100.0)
        mirna_target_interactions.rename(columns={"weighted context++ score percentile": "weight",
                                                  "MiRBase ID": "MIR",
                                                  "Gene Symbol": "GE"}, inplace=True)
        mir_target_network = nx.from_pandas_edgelist(mirna_target_interactions, source="MIR", target="GE",
                                                     edge_attr="weight", create_using=nx.DiGraph())
        return mir_target_network.edges(data=True)


    def process_genes_info(self):
        self.gene_info = pd.DataFrame(index=self.get_genes_list())

        self.gene_info.index.name = "MiRBase ID"
        self.gene_info = self.gene_info.join(self.targetScan_family_df.groupby("MiRBase ID").first(), on="MiRBase ID",how="left")

        self.gene_info = self.gene_info.join(self.HUGO_miRNA_gene_info_df, on="MiRBase ID")

        self.gene_info["Disease association"] = self.gene_info.index.map(
            self.mirnadisease.groupby("miRNA name")["Disease name"].apply('|'.join).to_dict())

        self.gene_info["locus_type"] = "RNA, micro"

        self.gene_info.rename(columns={'Mature sequence': 'Transcript sequence'}, inplace=True)
        self.gene_info["Transcript length"] = self.gene_info["Transcript sequence"].apply(
            lambda x: len(x) if type(x) is str else None)

    def get_genes_info(self):
        return self.gene_info


class CopyNumberVariation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "copyNumber.txt")
        super().__init__(cancer_type, file_path)


class DNAMethylation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "methylation_450.txt")
        super().__init__(cancer_type, file_path)


class ProteinExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, log2_transform=True):
        file_path = os.path.join(folder_path, "protein_RPPA.txt")
        super().__init__(cancer_type, file_path, log2_transform=log2_transform)

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_edgelist(HPRD_PPI, source=0, target=3,
                                          create_using=nx.DiGraph())

    def get_HPRD_PPI_network_edgelist(self):
        return self.HPRD_PPI_network.edges()

    def process_STRING_PPI_network(self, ppi_data_file_path):
        pass


