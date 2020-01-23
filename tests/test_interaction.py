from openomics.database import LncBase, MiRTarBase, STRING
from .test_multiomics import *

@pytest.fixture
def generate_LncBase():
    return LncBase("/data/datasets/Bioinformatics_ExternalData/lncBase/", )


@pytest.fixture
def generate_MiRTarBase():
    return MiRTarBase("/data/datasets/Bioinformatics_ExternalData/miRTarBase/", target_index="Target Gene")


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
