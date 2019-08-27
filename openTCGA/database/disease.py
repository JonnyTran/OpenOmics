
from openTCGA.database.annotation import *

class DiseaseAssociation(Database):

    @abstractmethod
    def get_disease_assocs(self, index): raise NotImplementedError
