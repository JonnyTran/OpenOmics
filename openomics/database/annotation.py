import os
from typing import List, Dict, Union
from abc import ABCMeta, abstractmethod
from io import StringIO
from os.path import expanduser
import difflib

import pandas
import pandas as pd
from Bio import SeqIO
from bioservices import BioMart

import openomics
from openomics.utils import GTF
from openomics.utils.df import concat_uniques_agg
from openomics.utils.io import mkdirs

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openomics")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), ".openomics", "databases")


class Dataset:
    COLUMNS_RENAME_DICT = None  # Needs initialization since subclasses may use this field

    def __init__(self, import_folder, file_resources=None, col_rename=None, **kwargs):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotations, functional annotations, sequences, and disease associations.
        Args:
            import_folder (str):
                The folder path containing the data files
            file_resources (dict): default None,
                Used to list required files for preprocessing of the database. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            col_rename (dict): default None,
                A dictionary to rename columns in the data table. If None, then automatically load defaults.
            **kwargs: Additional arguments that may be passed to load_data function
        """
        if not os.path.isdir(import_folder) or not os.path.exists(import_folder):
            raise IOError(import_folder)
        else:
            for _, filepath in file_resources.items():
                if not os.path.exists(filepath):
                    raise IOError(filepath)

        self.import_folder = import_folder
        self.file_resources = file_resources
        self.df = self.load_dataframe(file_resources, **kwargs)
        if col_rename is not None:
            self.df.rename(columns=col_rename, inplace=True)
        print("{}: {}".format(self.name(), self.df.columns.tolist()))

    @classmethod
    def name(cls):
        return cls.__name__

    def list_databases(self):
        return DEFAULT_LIBRARIES

    def get_annotations(self, index, columns):
        # type: (str, List[str]) -> pd.DataFrame
        """
        Returns the Database's DataFrame such that it's indexed by :param index:, which then applies a groupby operation
        and aggregates all other columns by concatenating all unique values.

        operation aggregates
        Args:
            index (str): The index column name of the Dataframe
            columns (list): a list of column names

        Returns:
            df (DataFrame): A dataframe to be used for annotation

        """
        if columns is not None:
            if index in columns:
                df = self.df.filter(items=columns)
                columns.pop(columns.index(index))
            else:
                df = self.df.filter(items=columns + [index])
        else:
            raise Exception("The columns argument must be a list such that it's subset of the following columns in the dataframe",
                            self.df.columns.tolist())

        if index != self.df.index.name and index in self.df.columns:
            df.set_index(index, inplace=True)

        # Groupby index, and Aggregate by all columns by concatenating unique values
        df = df.groupby(index).agg({k:concat_uniques_agg for k in columns})

        if df.index.duplicated().sum() > 0:
            raise ValueError("DataFrame must not have duplicates in index")
        return df

    @abstractmethod
    def load_dataframe(self, file_resources, **kwargs):
        # type: (dict, **str) -> pd.DataFrame
        """
        Handles data preprocessing given the file_resources input, and returns a DataFrame.

        Args:
            file_resources (dict): A dict with keys as filenames and values as full file path.
            **kwargs: Optional
        """
        raise NotImplementedError

    @abstractmethod
    def get_rename_dict(self, from_index, to_index):
        """
        Used to retrieve a lookup dictionary to convert from one index to another, e.g., gene_id to gene_name

        Args:
            from_index: an index on the DataFrame for key
            to_index: an index on the DataFrame for value

        Returns
            rename_dict (dict): a rename dict
        """
        raise NotImplementedError

    @abstractmethod
    def get_sequences(self, index, omic=None):
        """
        Returns a dictionary where keys are
        Args:
            omic (str): {"lncRNA", "microRNA", "messengerRNA"}
            index (str): {"gene_id", "gene_name", "transcript_id", "transcript_name"}
                The index
        """
        raise NotImplementedError



class Annotatable:
    """
    This abstract class provides an interface for the omics to annotate external data downloaded from various databases. These data will be imported as attribute information to the genes, or interactions between the genes.
    """
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    def get_annotations(self):
        if hasattr(self, "annotations"):
            return self.annotations
        else:
            raise Exception("Must run initialize_annotations() first.")

    def initialize_annotations(self, gene_list, index):
        if gene_list is None:
            gene_list = self.get_genes_list()

        self.annotations = pd.DataFrame(index=gene_list)
        self.annotations.index.name = index

    def annotate_genomics(self, database, index, columns, fuzzy_match=False):
        # type: (Dataset, str, List[str], bool) -> None
        """
        Performs a left outer join between the annotations and Database's DataFrame, on the index key. The index argument must be column present in both DataFrames.
        If there exists overlapping column in the join, then the fillna() is used to fill NaN values in the old column with non-NaN values from the new column.
        Args:
            database (openomics.annotation.Database): Database which contains annotations
            index (str): The column name which exists in both the annotations and Database's DataFrame
            columns (list): a list of column name to join to the annotations
            fuzzy_match (bool): default False. Whether to join the annotations by applying a fuzzy match on the index with difflib.get_close_matches(). It is very computationally expensive and thus should only be used sparingly.
        """
        database_annotations = database.get_annotations(index, columns)
        if fuzzy_match:
            database_annotations.index = database_annotations.index.map(lambda x: difflib.get_close_matches(x, self.annotations.index)[0])

        if index == self.annotations.index.name:
            self.annotations = self.annotations.join(database_annotations, on=index, rsuffix="_")
        else:
            old_index = self.annotations.index.name
            self.annotations = self.annotations.reset_index()
            self.annotations.set_index(index, inplace=True)
            self.annotations = self.annotations.join(database_annotations, on=index, rsuffix="_")
            self.annotations = self.annotations.reset_index()
            self.annotations.set_index(old_index, inplace=True)

        # Merge columns if the database DataFrame has overlapping columns with existing column
        duplicate_columns = [col for col in self.annotations.columns if col[-1] == "_"]
        for col in duplicate_columns:
            self.annotations[col.strip("_")].fillna(self.annotations[col], inplace=True, axis=0)
            self.annotations.drop(columns=col, inplace=True)

    def annotate_sequences(self, database, index, omic):
        # type: (Dataset, str, str) -> None
        self.annotations["Transcript sequence"] = self.annotations.index.map(
            database.get_sequences(index=index, omic=omic))

    def annotate_interactions(self, database, index):
        # type: (Dataset, str) -> None
        raise NotImplementedError

    def annotate_diseases(self, database, index):
        # type: (Dataset, str) -> None
        raise NotImplementedError


class RNAcentral(Dataset):
    COLUMNS_RENAME_DICT = {'ensembl_gene_id': 'gene_id',
                           'gene symbol': 'gene_name',
                           'external id': 'transcript_id',
                           'GO terms': 'go_id'}

    def __init__(self, import_folder, file_resources=None, col_rename=None, species=9606):
        self.species = species

        if file_resources is None:
            file_resources = {}
            file_resources["rnacentral_rfam_annotations.tsv"] = os.path.join(import_folder,
                                                                          "rnacentral_rfam_annotations.tsv")
            file_resources["gencode.tsv"] = os.path.join(import_folder, "gencode.tsv")

        if col_rename is None:
            col_rename = self.COLUMNS_RENAME_DICT

        super(RNAcentral, self).__init__(import_folder, file_resources, col_rename)

    def load_dataframe(self, file_resources, **kwargs):
        go_terms = pd.read_table(file_resources["rnacentral_rfam_annotations.tsv"],
                                 low_memory=True, header=None, names=["RNAcentral id", "GO terms", "Rfams"])
        go_terms["RNAcentral id"] = go_terms["RNAcentral id"].str.split("_", expand=True)[0]

        gencode_id = pd.read_table(file_resources["gencode.tsv"],
                                   low_memory=True, header=None,
                                   names=["RNAcentral id", "database", "external id", "species", "RNA type",
                                          "gene symbol"])

        gencode_id["species"] = gencode_id["species"].astype("O")
        if self.species is not None:
            gencode_id = gencode_id[gencode_id["species"] == self.species]

        lnc_go_terms = go_terms[go_terms["RNAcentral id"].isin(gencode_id["RNAcentral id"])].groupby("RNAcentral id")[
            "GO terms"].apply(lambda x: "|".join(x.unique()))
        lnc_rfams = go_terms[go_terms["RNAcentral id"].isin(gencode_id["RNAcentral id"])].groupby("RNAcentral id")[
            "Rfams"].apply(lambda x: "|".join(x.unique()))

        gencode_id["GO terms"] = gencode_id["RNAcentral id"].map(lnc_go_terms.to_dict())
        gencode_id["Rfams"] = gencode_id["RNAcentral id"].map(lnc_rfams.to_dict())
        gencode_id = gencode_id[gencode_id["GO terms"].notnull() | gencode_id["Rfams"].notnull()]

        return gencode_id


class GENCODE(Dataset):
    def __init__(self, import_folder, file_resources=None, col_rename=None, import_sequences="all",
                 replace_U2T=True):
        if file_resources is None:
            file_resources = {}
            file_resources["long_noncoding_RNAs.gtf"] = os.path.join(import_folder, "gencode.v29.long_noncoding_RNAs.gtf")
            file_resources["lncRNA_transcripts.fa"] = os.path.join(import_folder, "gencode.v29.lncRNA_transcripts.fa")
            file_resources["transcripts.fa"] = os.path.join(import_folder, "gencode.v29.transcripts.fa")

        self.import_sequences = import_sequences
        self.replace_U2T = replace_U2T

        super(GENCODE, self).__init__(import_folder, file_resources, col_rename=col_rename)

    def load_dataframe(self, file_resources, **kwargs):
        # Parse lncRNA gtf
        df = GTF.dataframe(file_resources["long_noncoding_RNAs.gtf"])
        df['gene_id'] = df['gene_id'].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
        df['transcript_id'] = df['transcript_id'].str.replace("[.].*", "")
        return df

    def get_sequences(self, index, omic=None):
        # Parse lncRNA & mRNA fasta
        if omic == openomics.transcriptomics.MessengerRNA.name():
            fasta_file = self.file_resources["transcripts.fa"]
        elif omic == openomics.transcriptomics.LncRNA.name():
            fasta_file = self.file_resources["lncRNA_transcripts.fa"]
        else:
            raise Exception("The omic argument must be one of the omic names")

        seq_dict = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            if index == "gene_id":
                key = record.id.split("|")[1].split(".")[0] # gene id
            elif index == "gene_name":
                key = record.id.split("|")[5]  # gene_name
            elif index == "transcript_id":
                key = record.id.split("|")[0].split(".")[0]  # transcript ID
            elif index == "transcript_name":
                key = record.id.split("|")[4]  # transcript_name
            else:
                raise Exception("The level argument must be one of 'gene_id', 'transcript_id', or 'gene_name', or 'transcript_name'")

            sequence_str = str(record.seq)
            if self.replace_U2T:
                sequence_str = sequence_str.replace("U", "T")

            # If index by gene, then select transcript sequences either by "shortest", "longest" or "all"
            if "gene" in index and self.import_sequences == "shortest":
                if key not in seq_dict:
                    seq_dict[key] = sequence_str
                else:
                    if len(seq_dict[key]) > len(sequence_str):
                        seq_dict[key] = sequence_str
            elif "gene" in index and self.import_sequences == "longest":
                if key not in seq_dict:
                    seq_dict[key] = sequence_str
                else:
                    if len(seq_dict[key]) < len(sequence_str):
                        seq_dict[key] = sequence_str
            elif "gene" in index and self.import_sequences == "all":
                if key not in seq_dict:
                    seq_dict[key] = [sequence_str, ]
                else:
                    seq_dict[key].append(sequence_str)
            else:
                seq_dict[key] = sequence_str

        return seq_dict

    def get_rename_dict(self, from_index, to_index):
        ensembl_id_to_gene_name = pd.Series(self.df['gene_name'].values,
                                            index=self.df['gene_id']).to_dict()
        return ensembl_id_to_gene_name


class MirBase(Dataset):
    def __init__(self, import_folder, RNAcentral_folder, file_resources=None, col_rename=None,
                 species=9606, import_sequences="all", replace_U2T=True):
        """

        Args:
            import_folder:
            RNAcentral_folder:
            file_resources:
            col_rename:
            species:
            import_sequences (str): {"longest", "shortest", "all"}
                Whether to select the longest, shortest, or a list of all transcript sequences when aggregating transcript sequences by gene_id or gene_name.
            replace_U2T:
        """
        if file_resources is None:
            file_resources = {}
            file_resources["aliases.txt"] = os.path.join(import_folder, "aliases.txt")
            file_resources["mature.fa"] = os.path.join(import_folder, "mature.fa")
            file_resources["rnacentral.mirbase.tsv"] = os.path.join(RNAcentral_folder, "mirbase.tsv")
            file_resources["rnacentral_rfam_annotations.tsv"] = os.path.join(RNAcentral_folder, "rnacentral_rfam_annotations.tsv")

        self.import_sequences = import_sequences
        self.replace_U2T = replace_U2T
        self.species = species
        super(MirBase, self).__init__(import_folder, file_resources, col_rename)

    def load_dataframe(self, file_resources, **kwargs):
        rnacentral_mirbase = pd.read_table(file_resources["rnacentral.mirbase.tsv"], low_memory=True, header=None,
                                   names=["RNAcentral id", "database", "mirbase id", "species", "RNA type", "gene name"],
                                   # dtype="O",
                                   index_col="mirbase id")
        #
        rnacentral_mirbase["species"] = rnacentral_mirbase["species"].astype("O")
        if self.species is not None:
            rnacentral_mirbase = rnacentral_mirbase[rnacentral_mirbase["species"] == self.species]

        mirbase_aliases = pd.read_table(file_resources["aliases.txt"], low_memory=True, header=None,
                                     names=["mirbase id", "gene_name"], dtype="O")
        mirbase_aliases = mirbase_aliases.join(rnacentral_mirbase, on="mirbase id", how="inner")

        # # Expanding miRNA names in each MirBase Ascension ID
        mirna_names = mirbase_aliases.apply(lambda x: pd.Series(x['gene_name'].split(";")[:-1]), axis=1).stack().reset_index(
            level=1, drop=True)
        mirna_names.name = "gene_name"
        mirbase_aliases = mirbase_aliases.drop('gene_name', axis=1).join(mirna_names)

        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.lower()
        # mirbase_name["miRNA name"] = mirbase_name["miRNA name"].str.replace("-3p.*|-5p.*", "")

        return mirbase_aliases

    def get_sequences(self, index="gene_name", omic=None):
        seq_dict = {}

        for record in SeqIO.parse(self.file_resources["mature.fa"], "fasta"):
            gene_name = str(record.id)

            sequence_str = str(record.seq)
            if self.replace_U2T:
                sequence_str = sequence_str.replace("U", "T")

            if self.import_sequences == "shortest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) > len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif self.import_sequences == "longest":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = sequence_str
                else:
                    if len(seq_dict[gene_name]) < len(sequence_str):
                        seq_dict[gene_name] = sequence_str
            elif self.import_sequences == "all":
                if gene_name not in seq_dict:
                    seq_dict[gene_name] = [sequence_str, ]
                else:
                    seq_dict[gene_name].append(sequence_str)
            else:
                seq_dict[gene_name] = sequence_str

        return seq_dict


class BioMartManager:
    __metaclass__ = ABCMeta

    def __init__(self, dataset, attributes, host, filename):
        pass

    def query_biomart(self, dataset, attributes, host="www.ensembl.org", cache=True, save_filename=None):
        bm = BioMart(host=host)
        bm.new_query()
        bm.add_dataset_to_xml(dataset)
        for at in attributes:
            bm.add_attribute_to_xml(at)
        xml_query = bm.get_xml()

        print("Querying {} from {} with attributes {}...".format(dataset, host, attributes))
        results = bm.query(xml_query)
        df = pd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", index_col=None, low_memory=True)

        if cache:
            self.cache_dataset(dataset, df, save_filename)
        return df

    def cache_dataset(self, dataset, dataframe, save_filename):
        if save_filename is None:
            mkdirs(DEFAULT_CACHE_PATH)
            save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(dataset))
        dataframe.to_csv(save_filename, sep="\t", index=False)
        return save_filename

    def retrieve_dataset(self, host, dataset, attributes, filename):
        filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(filename))
        if os.path.exists(filename):
            df = pd.read_csv(filename, sep="\t", low_memory=True)
        else:
            df = self.query_biomart(host=host, dataset=dataset, attributes=attributes,
                                    cache=True, save_filename=filename)
        return df


class EnsemblGenes(BioMartManager, Dataset):
    COLUMNS_RENAME_DICT = {'ensembl_gene_id': 'gene_id',
                           'external_gene_name': 'gene_name',
                           'ensembl_transcript_id': 'transcript_id',
                           'external_transcript_name': 'transcript_name',
                           'rfam': 'Rfams'}

    def __init__(self, dataset="hsapiens_gene_ensembl",
                 attributes=None,
                 host="www.ensemble.org", filename=False):
        if attributes is None:
            attributes = ['ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id',
                          'external_transcript_name',
                          'chromosome_name', 'transcript_start', 'transcript_end', 'transcript_length',
                          'gene_biotype', 'transcript_biotype',
                          'rfam', 'go_id', ]
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)

        self.df.rename(columns=self.COLUMNS_RENAME_DICT,
                       inplace=True)
        print(self.name(), self.df.columns.tolist())

    def load_dataframe(self, datasets, attributes, host, filename=None, ):
        return self.retrieve_dataset(host, datasets, attributes, filename)

    def get_rename_dict(self, from_index="gene_id", to_index="gene_name"):
        geneid_to_genename = self.df[self.df[to_index].notnull()]\
            .groupby(from_index)[to_index]\
            .apply(concat_uniques_agg).to_dict()
        return geneid_to_genename

    def get_functional_annotations(self, omic, index):
        geneid_to_go = self.df[self.df["go_id"].notnull()]\
            .groupby(index)["go_id"]\
            .apply(lambda x: "|".join(x.unique())).to_dict()
        return geneid_to_go


class EnsemblGeneSequences(EnsemblGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl",
                 attributes=None,
                 host="www.ensemble.org", filename=False):
        if attributes is None:
            attributes = ['ensembl_gene_id', 'gene_exon_intron', 'gene_flank', 'coding_gene_flank', 'gene_exon',
                          'coding']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, filename=self.filename, host=self.host,
                                      attributes=attributes, )
        self.df.rename(columns=self.COLUMNS_RENAME_DICT,
                       inplace=True)


class EnsemblTranscriptSequences(EnsemblGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl",
                 attributes=None,
                 host="www.ensemble.org", filename=False):
        if attributes is None:
            attributes = ['ensembl_transcript_id', 'transcript_exon_intron', 'transcript_flank',
                          'coding_transcript_flank',
                          '5utr', '3utr']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)
        self.df.rename(columns=self.COLUMNS_RENAME_DICT,
                       inplace=True)


class EnsemblSNP(EnsemblGenes):
    def __init__(self, dataset="hsapiens_snp",
                 attributes=None,
                 host="www.ensemble.org", filename=False):
        if attributes is None:
            attributes = ['synonym_name', 'variation_names', 'minor_allele',
                          'associated_variant_risk_allele',
                          'ensembl_gene_stable_id', 'ensembl_transcript_stable_id',
                          'phenotype_name',
                          'chr_name', 'chrom_start', 'chrom_end']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)


class EnsemblSomaticVariation(EnsemblGenes):
    def __init__(self, dataset="hsapiens_snp_som",
                 attributes=None,
                 host="www.ensemble.org", filename=False):
        if attributes is None:
            attributes = ['somatic_variation_name', 'somatic_source_name', 'somatic_allele', 'somatic_minor_allele',
                          'somatic_clinical_significance', 'somatic_validated', 'somatic_transcript_location',
                          'somatic_mapweight',
                          'somatic_chromosome_start', 'somatic_chromosome_end']
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.host = host
        self.df = self.load_dataframe(datasets=dataset, attributes=attributes, host=self.host,
                                      filename=self.filename)


class NONCODE(Dataset):
    # TODO need more fix
    def __init__(self, import_folder, file_resources=None, col_rename=None, **kwargs):
        if file_resources is None:
            file_resources = {}
            file_resources["NONCODEv5_source"] = os.path.join(import_folder, "NONCODEv5_source")
            file_resources["NONCODEv5_Transcript2Gene"] = os.path.join(import_folder, "NONCODEv5_Transcript2Gene")
            file_resources["NONCODEv5_human.func"] = os.path.join(import_folder, "NONCODEv5_human.func")

        super().__init__(import_folder, file_resources, col_rename, **kwargs)

    def load_dataframe(self, file_resources, **kwargs):
        source_df = pd.read_table(file_resources["NONCODEv5_source"], header=None)
        source_df.columns = ["NONCODE Transcript ID", "name type", "Gene ID"]

        transcript2gene_df = pd.read_table(file_resources["NONCODEv5_Transcript2Gene"], header=None)
        transcript2gene_df.columns = ["NONCODE Transcript ID", "NONCODE Gene ID"]

        self.noncode_func_df = pd.read_table(file_resources["NONCODEv5_human.func"], header=None)
        self.noncode_func_df.columns = ["NONCODE Gene ID", "GO terms"]
        self.noncode_func_df.set_index("NONCODE Gene ID", inplace=True)

        # Convert to NONCODE transcript ID for the functional annotation data
        self.noncode_func_df["NONCODE Transcript ID"] = self.noncode_func_df.index.map(
            pd.Series(transcript2gene_df['NONCODE Transcript ID'].values,
                      index=transcript2gene_df['NONCODE Gene ID']).to_dict())

        # Convert NONCODE transcript ID to gene names
        source_gene_names_df = source_df[source_df["name type"] == "NAME"].copy()

        self.noncode_func_df["Gene Name"] = self.noncode_func_df["NONCODE Transcript ID"].map(
            pd.Series(source_gene_names_df['Gene ID'].values,
                      index=source_gene_names_df['NONCODE Transcript ID']).to_dict())


# Constants
DEFAULT_LIBRARIES=["10KImmunomes"
"BioGRID"
"CCLE"
"DisGeNET"
"ENSEMBL"
"GENCODE"
"GeneMania"
"GeneOntology"
"GlobalBiobankEngine"
"GTEx"
"HMDD_miRNAdisease"
"HPRD_PPI"
"HUGO_Gene_names"
"HumanBodyMapLincRNAs"
"IntAct"
"lncBase"
"LNCipedia"
"LncReg"
"lncRInter"
"lncrna2target"
"lncRNA_data_repository"
"lncrnadisease"
"lncRNome"
"mirbase"
"miRTarBase"
"NHLBI_Exome_Sequencing_Project"
"NONCODE"
"NPInter"
"PIRD"
"RegNetwork"
"RISE_RNA_Interactions"
"RNAcentral"
"StarBase_v2.0"
"STRING_PPI"
"TargetScan"]
