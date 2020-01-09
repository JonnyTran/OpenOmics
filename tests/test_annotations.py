import pytest

from openomics.database import GENCODE, RNAcentral, MirBase


@pytest.fixture
def generate_GENCODE_ftp():
    return GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                   file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                   "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                   "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                   import_sequences="shortest")


@pytest.fixture
def generate_RNACentral_ftp():
    return RNAcentral(path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/",
                      file_resources={
                          "rnacentral_rfam_annotations.tsv": "go_annotations/rnacentral_rfam_annotations.tsv.gz",
                          "gencode.tsv": "id_mapping/database_mappings/gencode.tsv"},
                      )


@pytest.fixture
def generate_MirBase_ftp():
    return MirBase(path="ftp://mirbase.org/pub/mirbase/CURRENT/")


# @pytest.fixture
# def generate_GTEx_expressions():
#     return GTEx(path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/")


def test_import_gencode_db(generate_GENCODE_ftp):
    assert generate_GENCODE_ftp.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/'


def test_import_rnacentral_db(generate_RNACentral_ftp):
    assert generate_RNACentral_ftp.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/'


def test_import_mirbase_db(generate_MirBase_ftp):
    assert generate_MirBase_ftp.data_path == "ftp://mirbase.org/pub/mirbase/CURRENT/"

# def test_import_GTEx(generate_GTEx_expressions):
#     assert generate_GTEx_expressions.data_path == "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
