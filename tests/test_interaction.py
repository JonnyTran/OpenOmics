import pytest

from openomics.database import LncBase, MiRTarBase


@pytest.fixture
def generate_LncBase():
    return LncBase("/data/datasets/Bioinformatics_ExternalData/lncBase/", )


@pytest.fixture
def generate_MiRTarBase():
    return MiRTarBase("/data/datasets/Bioinformatics_ExternalData/miRTarBase/", target_index="Target Gene")
