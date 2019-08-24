import os
from abc import ABCMeta, abstractmethod
from io import StringIO
from os.path import expanduser

import pandas as pd
from Bio import SeqIO
from bioservices import BioMart

from openTCGA.utils import GTF
from openTCGA.utils.io import mkdirs
from openTCGA.utils.df import concat_uniques_agg

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openTCGA")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), ".openTCGA", "databases")


class Database:
    __metaclass__ = ABCMeta
    def query_biomart(self, dataset, attributes, host="www.ensembl.org", cache=True, save_filename=None):
        bm = BioMart(host=host)
        bm.new_query()
        bm.add_dataset_to_xml(dataset)
        for at in attributes:
            bm.add_attribute_to_xml(at)
        xml_query = bm.get_xml()

        print("Querying {} from {} with attributes {}...".format(dataset, host, attributes))
        results = bm.query(xml_query)
        df = pd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", index_col=None)

        if cache:
            self.cache_dataset(dataset, df, save_filename)
        return df

    def cache_dataset(self, dataset, dataframe, save_filename):
        if save_filename is None:
            mkdirs(DEFAULT_CACHE_PATH)
            save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(dataset))
        dataframe.to_csv(save_filename, sep="\t", index=False)
        return save_filename

    def retrieve_dataset(self, dataset, attributes, filename):
        filename = os.path.join(DEFAULT_CACHE_PATH, "{}.tsv".format(filename))
        if os.path.exists(filename):
            df = pd.read_csv(filename, sep="\t", low_memory=False)
        else:
            df = self.query_biomart(host="www.ensembl.org", dataset=dataset, attributes=attributes,
                                    cache=True, save_filename=filename)
        return df

    @classmethod
    def list_databases(cls,):
        return DEFAULT_LIBRARIES

    @abstractmethod
    def load_datasets(self, datasets, filename, **args):
        raise NotImplementedError

    @abstractmethod
    def get_id2name_dict(self, from_index, to_index) -> dict: raise NotImplementedError
    @abstractmethod
    def get_genomic_annotations(self, modality, index, columns): raise NotImplementedError
    @abstractmethod
    def get_functional_annotations(self, modality, index): raise NotImplementedError
    @abstractmethod
    def get_sequences(self, modality, index, *arg) -> dict: raise NotImplementedError
    @abstractmethod
    def get_interactions(self, modality, index): raise NotImplementedError
    @abstractmethod
    def get_disease_assocs(self, index): raise NotImplementedError

class Annotatable:
    __metaclass__ = ABCMeta

    def get_annotations(self):
        if hasattr(self, "annotations"):
            return self.annotations
        else:
            raise Exception("Must run initialize_annotations() first.")


    def get_network_edgelist(self):
        if hasattr(self, "network"):
            return self.network.edges(data=True)
        else:
            print(self.__class__.__str__(), "does not have network interaction data yet. (at self.network)")
            return None

    @abstractmethod
    def initialize_annotations(self, gene_list=None, index=None): raise NotImplementedError

    @abstractmethod
    def annotate_genomics(self, database: Database, index, columns):
        """
        Returns a DataFrame indexed by :index:. There must be no duplicates in the index column.

        :param database:
        :param index:
        :param columns:
        """
        raise NotImplementedError

    @abstractmethod
    def annotate_functions(self, database: Database, index, columns): raise NotImplementedError

    @abstractmethod
    def annotate_sequences(self, database: Database, index, **kwargs): raise NotImplementedError

    @abstractmethod
    def annotate_interactions(self, database: Database, index): raise NotImplementedError

    @abstractmethod
    def annotate_diseases(self, database: Database, index): raise NotImplementedError


class GENCODE(Database):
    def __init__(self, import_folder=None, version="v29", modalities=["GE", "LNC"], import_sequences="shortest", replace_U2T=True) -> None:
        if import_folder is not None:
            if not os.path.isdir(import_folder) or not os.path.exists(import_folder):
                raise NotADirectoryError(import_folder)
            self.folder_path = import_folder
        else:
            self.folder_path = os.path.join(DEFAULT_LIBRARY_PATH, "GENCODE")

            if os.path.exists(self.folder_path):
                print("Default library path + GENCODE doesn't exist, importing GENCODE datasets")
                # TODO import GENCODE database
                print(os.listdir(self.folder_path))

        self.file_resources = {}
        self.file_resources["long_noncoding_RNAs.gtf"] = os.path.join(self.folder_path, "gencode.v29.long_noncoding_RNAs.gtf")
        self.file_resources["lncRNA_transcripts.fa"] = os.path.join(self.folder_path, "gencode.v29.lncRNA_transcripts.fa")
        self.file_resources["transcripts.fa"] = os.path.join(self.folder_path, "gencode.v29.transcripts.fa")

        for k, v in self.file_resources.items():
            if not os.path.exists(v):
                print(FileNotFoundError(v))

        self.import_sequences = import_sequences
        self.replace_U2T = replace_U2T
        self.load_datasets()

    def load_datasets(self, datasets=None, filename=None):
        # Parse lncRNA gtf
        self.GENCODE_LncRNA_info = GTF.dataframe(self.file_resources["long_noncoding_RNAs.gtf"])
        self.GENCODE_LncRNA_info['gene_id'] = self.GENCODE_LncRNA_info['gene_id'].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
        self.GENCODE_LncRNA_info['transcript_id'] = self.GENCODE_LncRNA_info['transcript_id'].str.replace("[.].*", "")

    def get_genomic_annotations(self, modality, index, columns):
        if modality == "LNC":
            df = self.GENCODE_LncRNA_info.set_index(index)
        else:
            raise NotImplementedError

        if columns is not None:
            df = df.filter(items=columns)

        return df

    def get_sequences(self, modality, index="gene_id"):
        # Parse lncRNA & mRNA fasta
        if modality == "GE":
            fasta_file = self.file_resources["transcripts.fa"]
        elif modality == "LNC":
            fasta_file = self.file_resources["lncRNA_transcripts.fa"]
        else:
            raise Exception("The modality argument must be one of 'LNC', 'GE'")

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

    def get_id2name_dict(self, left_index, right_index):
        if modality == "LNC":
            ensembl_id_to_gene_name = pd.Series(self.GENCODE_LncRNA_info['gene_name'].values,
                                                index=self.GENCODE_LncRNA_info['gene_id']).to_dict()
            return ensembl_id_to_gene_name


class EnsembleGenes(Database):
    def __init__(self, dataset="hsapiens_gene_ensembl", filename=None) -> None:
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.attributes = ['ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'external_transcript_name',
                      'rfam', 'go_id',
                      'chromosome_name', 'transcript_stwart', 'transcript_end', 'transcript_length',
                      'cds_start', 'cds_end', 'cds_length', '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end',
                      'gene_biotype', 'transcript_biotype']

        self.df = self.load_datasets(datasets=dataset, attributes=self.attributes, filename=self.filename)
        self.df.rename(columns={'ensembl_gene_id': 'gene_id',
                                'external_gene_name': 'gene_name',
                                'ensembl_transcript_id': 'transcript_id',
                                'external_transcript_name': 'transcript_name'},
                       inplace=True)

    def load_datasets(self, datasets, attributes, filename=None):
        return self.retrieve_dataset(datasets, attributes, filename)

    def get_id2name_dict(self, from_index="gene_id", to_index="gene_name"):
        geneid_to_genename = self.df[self.df[to_index].notnull()]\
            .groupby(from_index)[to_index]\
            .apply(concat_uniques_agg).to_dict()
        return geneid_to_genename

    def get_genomic_annotations(self, modality, index, columns):
        if columns is not None:
            df = self.df.filter(items=columns)
        else:
            df = self.df

        if index != self.df.index.name and index in self.df.columns:
            df.set_index(index, inplace=True)

        # Groupby index, and Aggregate by all columns by concatenating unique values
        df = df.groupby(index).agg({k:concat_uniques_agg for k in columns})

        if df.index.duplicated().sum() > 0:
            raise ValueError("DataFrame must not have duplicates in index")
        return df

    def get_functional_annotations(self, modality, index):
        geneid_to_go = self.df[self.df["go_id"].notnull()]\
            .groupby(index)["go_id"]\
            .apply(lambda x: "|".join(x.unique())).to_dict()
        return geneid_to_go


class EnsembleGeneSequences(EnsembleGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl") -> None:
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.attributes = ['ensembl_gene_id', 'gene_exon_intron', 'gene_flank', 'coding_gene_flank', 'gene_exon', 'coding']
        self.df = self.load_datasets(datasets=dataset, filename=self.filename, attributes=self.attributes, )
        self.df.rename(columns={'ensembl_gene_id': 'gene_id'},
                       inplace=True)
        
class EnsembleTranscriptSequences(EnsembleGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl", filename=None) -> None:
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.attributes = ['ensembl_transcript_id', 'transcript_exon_intron', 'transcript_flank', 'coding_transcript_flank',
                      '5utr', '3utr']
        self.df = self.load_datasets(datasets=dataset, attributes=self.attributes, filename=self.filename)
        self.df.rename(columns={'ensembl_transcript_id': 'transcript_id'},
                       inplace=True)

class EnsembleSNP(EnsembleGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl", filename=None) -> None:
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.attributes = ['variation_name', 'allele', 'minor_allele', 'mapweight', 'validated', 'allele_string_2076',
                      'clinical_significance',
                      'transcript_location', 'snp_chromosome_strand', 'chromosome_start', 'chromosome_end']
        self.df = self.load_datasets(datasets=dataset, attributes=self.attributes, filename=self.filename)

class EnsembleSomaticVariation(EnsembleGenes):
    def __init__(self, dataset="hsapiens_gene_ensembl", filename=None) -> None:
        self.filename = "{}.{}".format(dataset, self.__class__.__name__)
        self.attributes = ['somatic_variation_name', 'somatic_source_name', 'somatic_allele', 'somatic_minor_allele',
                      'somatic_clinical_significance', 'somatic_validated', 'somatic_transcript_location',
                      'somatic_mapweight',
                      'somatic_chromosome_start', 'somatic_chromosome_end']
        self.df = self.load_datasets(datasets=dataset, attributes=self.attributes, filename=self.filename)


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

