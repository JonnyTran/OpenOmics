import pytest

from openomics.database import GENCODE


@pytest.fixture
def generate_GENCODE_db_ftp():
    return GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                   file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                   "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                   "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                   import_sequences="shortest")


def test_import_expression_table_size(generate_TCGA_LUAD_MessengerRNA):
    assert generate_TCGA_LUAD_MessengerRNA.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/'
