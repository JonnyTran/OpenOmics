import os
import pandas as pd
import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from io import StringIO
from os.path import expanduser
from bioservices import BioMart
from pkg_resources import resource_filename

from openTCGA.utils.utils import mkdirs

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openTCGA")
DEFAULT_LIBRARY_PATH = os.path.join(expanduser("~"), "Bioinformatics_ExternalData")


class Annotator(object):
    def __init__(self, organism='human',) -> None:
        super().__init__()


def get_ensemble_genes(filename=None, dataset="hsapiens_gene_ensembl"):
    filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
    if os.path.exists(filename):
        df = pd.read_csv(filename, sep="\t")
    else:
        df = query_biomart(dataset=dataset)
    return df


def query_biomart(host="www.ensembl.org", dataset="hsapiens_gene_ensembl",
                  attributes=['ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'go_id'],
                  save_filename=None):
    bm = BioMart(host=host)

    # Start query
    bm.new_query()
    bm.add_dataset_to_xml(dataset)
    for at in attributes:
        bm.add_attribute_to_xml(at)
    xml_query = bm.get_xml()

    results = bm.query(xml_query)
    df = pd.read_csv(StringIO(results), header=None, names=attributes, sep="\t", index_col=None)

    if save_filename is None:
        mkdirs(DEFAULT_CACHE_PATH)
        save_filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
    df.to_csv(save_filename, sep="\t", index=False)
    return df


def retry(num=5):
    """"retry connection.

        define max tries num
        if the backoff_factor is 0.1, then sleep() will sleep for
        [0.1s, 0.2s, 0.4s, ...] between retries.
        It will also force a retry if the status code returned is 500, 502, 503 or 504.

    """
    s = requests.Session()
    retries = Retry(total=num, backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])
    s.mount('http://', HTTPAdapter(max_retries=retries))

    return s


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

