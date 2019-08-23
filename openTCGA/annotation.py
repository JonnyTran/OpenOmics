import os
from abc import ABCMeta, abstractmethod
import pandas as pd
from io import StringIO
from os.path import expanduser
from bioservices import BioMart
from Bio import SeqIO
from Bio.UniProt import GOA

from openTCGA.utils.io import mkdirs
from openTCGA.utils import GTF

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openTCGA")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), ".openTCGA", "databases")


class Database:
    __metaclass__ = ABCMeta
    @classmethod
    def query_biomart(cls, host="www.ensembl.org", dataset="hsapiens_gene_ensembl",
                      attributes=['ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'go_id'],
                      cache=True, save_filename=None):
        bm = BioMart(host=host)
        # Start query
        bm.new_query()
        bm.add_dataset_to_xml(dataset)
        for at in attributes:
            bm.add_attribute_to_xml(at)
        xml_query = bm.get_xml()

        print("Querying {} from {}...".format(dataset, host))
        results = bm.query(xml_query)
        df = pd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", index_col=None)

        if cache:
            cls.cache_database(dataset, df, save_filename)
        return df

    @classmethod
    def cache_database(cls, database, dataframe, save_filename):
        if save_filename is None:
            mkdirs(DEFAULT_CACHE_PATH)
            save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(database))
        dataframe.to_csv(save_filename, sep="\t", index=False)
        return save_filename

    @classmethod
    def retrieve_database(cls, dataset="hsapiens_gene_ensembl", filename=None):
        filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
        if os.path.exists(filename):
            df = pd.read_csv(filename, sep="\t")
        else:
            df = cls.query_biomart(dataset=dataset)

        return df

    @classmethod
    def list_databases(cls,):
        return DEFAULT_LIBRARIES

    @abstractmethod
    def load_datasets(self, datasets, filename): raise NotImplementedError

    @abstractmethod
    def genename(self) -> dict: raise NotImplementedError
    @abstractmethod
    def genomic_annotations(self, modality): raise NotImplementedError
    @abstractmethod
    def functional_annotations(self, modality): raise NotImplementedError
    @abstractmethod
    def sequences(self, modality, *arg) -> dict: raise NotImplementedError
    @abstractmethod
    def interactions(self, modality): raise NotImplementedError
    @abstractmethod
    def disease_associations(self): raise NotImplementedError

class Annotatable:
    __metaclass__ = ABCMeta
    # @classmethod
    # def version(self): return "1.0"
    @abstractmethod
    def annotate(self, database:Database, key, level): raise NotImplementedError


class GENCODE(Database):
    def __init__(self, import_folder=None, version="v29", sequence=["GE", "LNC", "MIR"], import_sequences="all", replace_U2T=True) -> None:
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

        # Prase lncRNA & mRNA fasta
        self.locus_type_dict = {}
        self.seq_dict = {}

        for modality in ["GE", "LNC"]:
            if modality == "GE":
                fasta_file = self.file_resources["transcripts.fa"]
            elif modality == "LNC":
                fasta_file = self.file_resources["lncRNA_transcripts.fa"]

            self.seq_dict[modality] = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                # gene_id = record.id.split("|")[1]
                gene_name = record.id.split("|")[5]

                sequence_str = str(record.seq)
                if self.replace_U2T:
                    sequence_str = sequence_str.replace("U", "T")
                if self.import_sequences == "shortest":
                    if gene_name not in self.seq_dict[modality]:
                        self.seq_dict[modality][gene_name] = sequence_str
                    else:
                        if len(self.seq_dict[modality][gene_name]) > len(sequence_str):
                            self.seq_dict[modality][gene_name] = sequence_str
                elif self.import_sequences == "longest":
                    if gene_name not in self.seq_dict[modality]:
                        self.seq_dict[modality][gene_name] = sequence_str
                    else:
                        if len(self.seq_dict[modality][gene_name]) < len(sequence_str):
                            self.seq_dict[modality][gene_name] = sequence_str
                elif self.import_sequences == "all":
                    if gene_name not in self.seq_dict[modality]:
                        self.seq_dict[modality][gene_name] = [sequence_str, ]
                    else:
                        self.seq_dict[modality][gene_name].append(sequence_str)
                else:
                    self.seq_dict[modality][gene_name] = sequence_str

                # add locus type for mRNAs
                if modality == self.file_resources["transcripts.fa"]:
                    if ~(gene_name in self.locus_type_dict):
                        self.locus_type_dict[gene_name] = record.id.split("|")[7]
                    else:
                        self.locus_type_dict[gene_name] = self.locus_type_dict[gene_name] + "|" + record.id.split("|")[7]

    def genomic_annotations(self, modality="LNC"):
        if modality == "LNC":
            return self.GENCODE_LncRNA_info
        elif modality == "GE":
            return self.locus_type_dict

    def sequences(self, modality):
        if modality == "LNC":
            return self.seq_dict[""]

    def genename(self):
        ensembl_id_to_gene_name = pd.Series(self.GENCODE_LncRNA_info['gene_name'].values,
                                            index=self.GENCODE_LncRNA_info['gene_id']).to_dict()
        return ensembl_id_to_gene_name

class Ensemble(Database):
    def __init__(self, import_folder=None) -> None:
        self.df = self.retrieve_database(dataset="hsapiens_gene_ensembl")

    def genename(self):
        geneid_to_genename = self.df[self.df["external_gene_name"].notnull()].groupby('ensembl_gene_id')["external_gene_name"].apply(lambda x: "|".join(x.unique())).to_dict()
        return geneid_to_genename

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

