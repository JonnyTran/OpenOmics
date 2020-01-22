import pytest

from openomics.database import LncBase, MiRTarBase, STRING


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
    df = generate_STRING.get_annotations(index="protein_name", columns=generate_STRING.df.columns)
    assert df.loc[0, "protein_name"] in generate_STRING.get_sequences(index="protein_name")
