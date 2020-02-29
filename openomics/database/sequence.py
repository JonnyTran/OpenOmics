from .base import Dataset
from abc import abstractmethod
from abc import abstractmethod

from .base import Dataset


class SequenceDataset(Dataset):
    def __init__(self, agg_sequences="all", replace_U2T=True, **kwargs):
        self.agg_sequences = agg_sequences
        self.replace_U2T = replace_U2T

        super(SequenceDataset, self).__init__(**kwargs)

    @abstractmethod
    def get_sequences(self, index, omic, agg_sequences, **kwargs):
        """
        Returns a dictionary where keys are
        Args:
            omic (str): {"lncRNA", "microRNA", "messengerRNA"}
            index (str): {"gene_id", "gene_name", "transcript_id", "transcript_name"}
                The index
        """
        raise NotImplementedError
