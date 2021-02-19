import dask.dataframe as dd

from openomics.database import GENCODE, MirBase
from .test_multiomics import *


@pytest.fixture
def generate_GENCODE():
    return GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                   file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                   "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz"}, )


@pytest.fixture
def generate_GENCODE_dask():
    return GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                   file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                   "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz"}, npartitions=8)


@pytest.fixture
def generate_MirBase_ftp():
    return MirBase(path="ftp://mirbase.org/pub/mirbase/CURRENT/")


def test_import_mirbase_db(generate_MirBase_ftp):
    """
    Args:
        generate_MirBase_ftp:
    """
    assert generate_MirBase_ftp.data_path == "ftp://mirbase.org/pub/mirbase/CURRENT/"


def test_import_GENCODE(generate_GENCODE):
    """
    Args:
        generate_GENCODE:
    """
    assert generate_GENCODE.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/'


def test_import_dask_GENCODE(generate_GENCODE_dask):
    assert isinstance(generate_GENCODE_dask.data, dd.DataFrame)


def test_annotate_GENCODE(generate_TCGA_LUAD, generate_GENCODE):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GENCODE:
    """
    generate_TCGA_LUAD.LncRNA.annotate_attributes(generate_GENCODE, on="gene_id", columns=['gene_name'])
    assert {'gene_name'}.issubset(
        generate_TCGA_LUAD.LncRNA.annotations.columns)


def test_annotate_dask_GENCODE(generate_TCGA_LUAD, generate_GENCODE_dask):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GENCODE_dask:
    """
    generate_GENCODE_dask.data = generate_GENCODE_dask.data[generate_GENCODE_dask.data["gene_id"].notnull()]

    generate_TCGA_LUAD.LncRNA.annotate_attributes(generate_GENCODE_dask,
                                                  on="gene_id",
                                                  columns=['gene_name', 'transcript_id'],
                                                  agg="concat")
    assert {'gene_name', 'transcript_id'}.issubset(
        generate_TCGA_LUAD.LncRNA.annotations.columns)


def test_annotate_sequence_GENCODE(generate_TCGA_LUAD, generate_GENCODE):
    generate_TCGA_LUAD.LncRNA.annotate_sequences(generate_GENCODE,
                                                 index="gene_id",
                                                 omic="LncRNA",
                                                 agg="longest")
    assert not generate_TCGA_LUAD.LncRNA.annotations["sequence"].empty
