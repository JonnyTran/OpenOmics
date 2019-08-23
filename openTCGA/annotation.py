import os
from abc import ABCMeta, abstractmethod
import pandas as pd
from io import StringIO
from os.path import expanduser
from bioservices import BioMart

from openTCGA.utils.utils import mkdirs

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openTCGA")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), "Bioinformatics_ExternalData")


class Database:
    __metaclass__ = ABCMeta
    @abstractmethod
    def load_annotations(self): raise NotImplementedError

class Annotatable:
    __metaclass__ = ABCMeta
    # @classmethod
    # def version(self): return "1.0"
    @abstractmethod
    def annotate(self, database:Database, key, level): raise NotImplementedError

class Annotation:
    def __init__(self, database:str, key:str, level:str) -> None:
        pass

class GeneAnnotation(Annotation):
    pass

class FunctionalAnnotation(Annotation):
    pass

class SequenceAnnotation(Annotation):
    pass

class DiseaseAssociation(Annotation):
    pass

class Interactions(Annotation):
    pass



def retrieve_database(dataset="hsapiens_gene_ensembl", filename=None):
    filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
    if os.path.exists(filename):
        ensemble_genes = pd.read_csv(filename, sep="\t")
    else:
        ensemble_genes = query_biomart(dataset=dataset)

    return ensemble_genes

def query_biomart(host="www.ensembl.org", dataset="hsapiens_gene_ensembl",
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
        cache_database(dataset, df, save_filename)
    return df


def cache_database(database, dataframe, save_filename):
    if save_filename is None:
        mkdirs(DEFAULT_CACHE_PATH)
        save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(database))
    dataframe.to_csv(save_filename, sep="\t", index=False)


# Constants
DEFAULT_LIBRARY=["10KImmunomes"
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

