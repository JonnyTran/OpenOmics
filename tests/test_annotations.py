import pytest

from openomics.database import GENCODE, RNAcentral


@pytest.fixture
def generate_GENCODE_db_ftp():
    return GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                   file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                   "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                   "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                   import_sequences="shortest")


@pytest.fixture
def generate_RNACentral_db_ftp():
    return RNAcentral(path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/",
                      file_resources={
                          "rnacentral_rfam_annotations.tsv": "go_annotations/rnacentral_rfam_annotations.tsv.gz",
                          "gencode.tsv": "id_mapping/database_mappings/gencode.tsv"},
                      )


def test_import_gencode_db(generate_GENCODE_db_ftp):
    assert generate_GENCODE_db_ftp.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/'


def test_import_rnacentral_db(generate_RNACentral_db_ftp):
    assert generate_RNACentral_db_ftp.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/'
