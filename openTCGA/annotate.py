import os
import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from pkg_resources import resource_filename
from os.path import expanduser

DEFAULT_CACHE_PATH = os.path.join(expanduser("~"), ".openTCGA")


class Annotator(object):
    def __init__(self, organism='human',) -> None:
        super().__init__()


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

