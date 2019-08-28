
from openomics.database.annotation import *


class DiseaseAssociation(Dataset):

    @abstractmethod
    def get_disease_assocs(self, index): raise NotImplementedError
