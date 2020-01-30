from openomics.database import MiRTarBase, STRING, LncRNA2Target
from .test_multiomics import *


@pytest.fixture
def generate_LncRNA2Target():
    return LncRNA2Target(version="low_throughput")


def test_import_LncRNA2Target(generate_LncRNA2Target):
    assert generate_LncRNA2Target.data_path == "http://123.59.132.21/lncrna2target/data/"


@pytest.fixture
def generate_MiRTarBase():
    return MiRTarBase()


def test_import_MiRTarBase(generate_MiRTarBase):
    assert generate_MiRTarBase.data_path == "http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/"


@pytest.fixture
def generate_STRING():
    return STRING(path="https://stringdb-static.org/download/",
                  species_id="9606",
                  source_col_name="protein1", target_col_name="protein2", )


def test_import_STRING(generate_STRING):
    assert generate_STRING.data_path == "https://stringdb-static.org/download/"


def test_annotate_STRING(generate_TCGA_LUAD, generate_STRING):
    generate_TCGA_LUAD.Protein.annotate_sequences(generate_STRING, index="protein_name")
    assert not generate_TCGA_LUAD.Protein.annotations["Transcript sequence"].empty
