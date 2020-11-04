from openomics.database.base import Dataset


class DiseaseAssociation(Dataset):

    @abstractmethod
    def get_disease_assocs(self, index): raise NotImplementedError
